/*============================================================================
 * Main functions dedicated to the modelling of two-phase flows in a porous
 * media. This media is always considered as unsaturated. Two sub-models are
 * considered: miscible (MTPF) or immiscible (ITPF)
 *============================================================================*/

/* VERS */

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
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_array.h"
#include "cs_cdovb_priv.h"
#include "cs_field.h"
#include "cs_gwf_priv.h"
#include "cs_gwf_soil.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_reco.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_tpf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf_tpf.c

  \brief Main functions dedicated to the modelling of two-phase flows in a
         porous media. This media is always considered as unsaturated. Two
         sub-models are considered: miscible (MTPF) or immiscible (ITPF).
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */

#define CS_GWF_TPF_DBG 0

/*============================================================================
 * Local definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling of miscible two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_mtpf_log_setup(cs_gwf_tpf_t   *mc)
{
  if (mc == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Water mass density: %5.3e, viscosity: %5.3e,"
                " molar mass: %5.3e\n", mc->l_mass_density, mc->l_viscosity,
                mc->w_molar_mass);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Henry constant: %5.3e\n", mc->henry_constant);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Gas component: viscosity: %5.3e, diffusivity"
                " in the liquid phase: %5.3e, molar mass: %5.3e\n",
                mc->g_viscosity, mc->l_diffusivity_h, mc->h_molar_mass);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the modelling of immiscible two-phase flows
 *
 * \param[in] mc   pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

static void
_itpf_log_setup(cs_gwf_tpf_t   *mc)
{
  if (mc == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Water mass density: %5.3e, viscosity: %5.3e\n",
                mc->l_mass_density, mc->l_viscosity);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Henry constant: %5.3e (Should be very low)\n",
                mc->henry_constant);
  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Gas component viscosity: %5.3e, molar mass: %5.3e\n",
                mc->g_viscosity, mc->h_molar_mass);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add a source term to the hg equation when a segregated approach
 *          is used. Case of a vertex-based scheme.
 *
 *          Generic function prototype for a hook during the cellwise building
 *          of the linear system
 *          Fit the cs_equation_build_hook_t prototype. This function may be
 *          called by different OpenMP threads
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context to cast for this discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] context     pointer to a context structure
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_build_hg_incr_st(const cs_equation_param_t     *eqp,
                  const cs_equation_builder_t   *eqb,
                  const void                    *eq_context,
                  const cs_cell_mesh_t          *cm,
                  void                          *context,
                  cs_hodge_t                    *mass_hodge,
                  cs_hodge_t                    *diff_hodge,
                  cs_cell_sys_t                 *csys,
                  cs_cell_builder_t             *cb)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(mass_hodge);

  const cs_cdovb_scaleq_t  *eqc = (const cs_cdovb_scaleq_t *)eq_context;
  assert(csys->n_dofs == cm->n_vc);

  cs_gwf_tpf_t  *mc = context;

  /* Modify the property data associated to the hodge operator related to the
     diffusion term */

  cs_property_data_t  *saved = diff_hodge->pty_data;
  cs_property_data_t  tmp = cs_property_data_define(saved->need_tensor,
                                                    saved->need_eigen,
                                                    mc->diff_hl_pty);

  diff_hodge->pty_data = &tmp;

  /* Diffusion part of the source term to add */

  cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag, diff_hodge);

  /* Define the local stiffness matrix: local matrix owned by the cellwise
     builder (store in cb->loc) */

  eqc->get_stiffness_matrix(cm, diff_hodge, cb);

  /* Retrieve the values pl^{k,n+1} */

  double  *vec = cb->values, *matvec = cb->values + cm->n_vc;
  for (int v = 0; v < cm->n_vc; v++)
    vec[v] = mc->l_pressure->val[cm->v_ids[v]];

  cs_sdm_square_matvec(cb->loc, vec, matvec);

  /* Update the rhs */

  for (int v = 0; v < cm->n_vc; v++)
    csys->rhs[v] -= matvec[v];

  /* Set the diffusion property data back to the initial pointer */

  diff_hodge->pty_data = saved;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        liquid phase when a two-phase flow model is used.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      pot_values   values to consider for the update
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cur2prev     true or false
 * \param[in, out] darcy        pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_l_darcy_update(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *cdoq,
                    const cs_real_t             *pot_values,
                    cs_real_t                    t_eval,
                    bool                         cur2prev,
                    cs_gwf_darcy_flux_t         *darcy)
{
  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;
  assert(adv != NULL);
  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = (cs_gwf_tpf_t *)darcy->update_input;
  cs_equation_t  *eq = mc->wl_eq;

  /* Update the array of flux values associated to the advection field.
   *
   * diff_pty of the wl_eq -> rho_l * abs_perm * krl / mu_l
   * Thus, one needs to divide by rho_l for the Darcy flux
   */

  cs_equation_compute_diffusive_flux(eq,
                                     NULL,
                                     pot_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  const double  inv_rho = 1./mc->l_mass_density;
  const cs_adjacency_t  *c2e = connect->c2e;
# pragma omp parallel if (6*cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2e->idx[cdoq->n_cells]; i++)
    darcy->flux_val[i] *= inv_rho;

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        gas phase when a two-phase flow model is used.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      pot_values   values to consider for the update
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cur2prev     true or false
 * \param[in, out] darcy        pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_g_darcy_update(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *cdoq,
                    const cs_real_t             *pot_values,
                    cs_real_t                    t_eval,
                    bool                         cur2prev,
                    cs_gwf_darcy_flux_t         *darcy)
{
  CS_UNUSED(connect);
  CS_UNUSED(cdoq);
  CS_UNUSED(cur2prev);

  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);
  assert(darcy->adv_field != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = darcy->update_input;
  cs_equation_t  *eq = mc->hg_eq;

  /* Update the array of flux values associated to the advection field */

  cs_equation_compute_diffusive_flux(eq,
                                     NULL,
                                     pot_values,
                                     darcy->flux_location,
                                     t_eval,
                                     darcy->flux_val);

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);

  /* Update the Darcy flux at boundary faces if needed */

  cs_field_t  *bdy_nflx =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

  if (bdy_nflx != NULL) { /* Values of the Darcy flux at boundary face exist */

    if (cur2prev)
      cs_field_current_to_previous(bdy_nflx);

    /* Set the new values of the field related to the normal boundary flux */

    cs_advection_field_across_boundary(adv, t_eval, bdy_nflx->val);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the advection field/arrays related to the Darcy flux in the
 *        gas phase when a two-phase flow model is used.
 *        This relies on the case of CDO-Vb schemes and a Darcy flux defined
 *        at dual faces
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      pot_values   values to consider for the update
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cur2prev     true or false
 * \param[in, out] darcy        pointer to the darcy flux structure
 */
/*----------------------------------------------------------------------------*/

static void
_tpf_t_darcy_update(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *cdoq,
                    const cs_real_t             *pot_values,
                    cs_real_t                    t_eval,
                    bool                         cur2prev,
                    cs_gwf_darcy_flux_t         *darcy)
{
  assert(darcy != NULL);
  assert(darcy->flux_val != NULL);
  assert(darcy->update_input != NULL);
  assert(darcy->adv_field != NULL);

  cs_adv_field_t  *adv = darcy->adv_field;

  if (adv->definition->type != CS_XDEF_BY_ARRAY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid definition of the advection field", __func__);

  cs_gwf_tpf_t  *mc = darcy->update_input;
  const cs_gwf_darcy_flux_t  *l_darcy = mc->l_darcy;
  const cs_gwf_darcy_flux_t  *g_darcy = mc->g_darcy;

  const double  hmh = mc->h_molar_mass * mc->henry_constant;
  const double  mh_ov_rt =
    mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);
  const cs_adjacency_t  *c2e = connect->c2e;

# pragma omp parallel if (6*cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2e->idx[cdoq->n_cells]; i++)
    darcy->flux_val[i] =
      hmh * l_darcy->flux_val[i] + mh_ov_rt * g_darcy->flux_val[i];

  /* Update the velocity field at cell centers induced by the Darcy flux */

  cs_field_t  *vel = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);
  assert(vel != NULL);
  if (cur2prev)
    cs_field_current_to_previous(vel);

  cs_advection_field_in_cells(adv, t_eval, vel->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new increment on the capillarity pressure
 *
 * \param[in]      mesh     pointer to a cs_mesh_t structure
 * \param[in]      wl_eq    pointer on tne "wl" equation structure
 * \param[in]      hg_eq    pointer on tne "wl" equation structure
 * \param[in, out] dpc_kp1  values of the increment on the capillarity pressure
 */
/*----------------------------------------------------------------------------*/

static void
_get_capillarity_pressure_increment(const cs_mesh_t       *mesh,
                                    const cs_equation_t   *wl_eq,
                                    const cs_equation_t   *hg_eq,
                                    cs_real_t             *dpc_kp1)
{
  const cs_equation_param_t  *wl_eqp = cs_equation_get_param(wl_eq);
  const cs_equation_param_t  *hg_eqp = cs_equation_get_param(hg_eq);
  const cs_equation_builder_t  *wl_eqb = cs_equation_get_builder(wl_eq);
  const cs_equation_builder_t  *hg_eqb = cs_equation_get_builder(hg_eq);
  const cs_real_t  *dpl_kp1 = wl_eqb->increment;
  const cs_real_t  *dpg_kp1 = hg_eqb->increment;

  if (wl_eqp->incremental_relax_factor < 1 ||
      hg_eqp->incremental_relax_factor < 1) {

    assert(wl_eqp->incremental_relax_factor > 0);
    assert(hg_eqp->incremental_relax_factor > 0);

    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      dpc_kp1[i] = hg_eqp->incremental_relax_factor*dpg_kp1[i]
                 - wl_eqp->incremental_relax_factor*dpl_kp1[i];

  }
  else {

    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++)
      dpc_kp1[i] = dpg_kp1[i] - dpl_kp1[i];

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, "%s: dpl=% 6.4e; dpg=% 6.4e\n", __func__,
                sqrt(cs_cdo_blas_square_norm_pvsp(dpl_kp1)),
                sqrt(cs_cdo_blas_square_norm_pvsp(dpg_kp1)));
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the convergence of the non-linear algorithm in case of TPF
 *         model. Work only with the capillarity pressure at vertices
 *
 * \param[in]       nl_algo_type type of non-linear algorithm
 * \param[in]       dpc_iter     cur. increment values for the capill. pressure
 * \param[in, out]  algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_nl_tpf_dpc_cvg(cs_param_nl_algo_t        nl_algo_type,
                      const cs_real_t          *dpc_iter,
                      cs_iter_algo_t           *algo)
{
  assert(algo != NULL);
  algo->prev_res = algo->res;

  double dpc_norm = cs_cdo_blas_square_norm_pvsp(dpc_iter);

  algo->res = dpc_norm;
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                possible divergence of the algorithm */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg_default(algo);

  if (algo->verbosity > 0) {

    if (algo->n_algo_iter == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                    cs_param_get_nl_algo_label(nl_algo_type));

    cs_log_printf(CS_LOG_DEFAULT,
                  "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                  cs_param_get_nl_algo_label(nl_algo_type),
                  algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the convergence of the non-linear algorithm in case of TPF
 *         model
 *
 * \param[in]       nl_algo_type type of non-linear algorithm
 * \param[in]       pg_pre_iter  previous iterate values for the gas pressure
 * \param[in, out]  pg_cur_iter  current iterate values for the gas pressure
 * \param[in]       pl_pre_iter  previous iterate values for the liquid pressure
 * \param[in, out]  pl_cur_iter  current iterate values for the liquid pressure
 * \param[in, out]  algo         pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_nl_tpf_cvg(cs_param_nl_algo_t        nl_algo_type,
                  const cs_real_t          *pg_pre_iter,
                  cs_real_t                *pg_cur_iter,
                  const cs_real_t          *pl_pre_iter,
                  cs_real_t                *pl_cur_iter,
                  cs_iter_algo_t           *algo)
{
  assert(algo != NULL);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON && algo->n_algo_iter > 0) {

    /* pg_* arrays gather pg and pl (this is done during the solve step) */

    cs_iter_algo_aa_update(algo,
                           pg_cur_iter, /* updated during the process */
                           pg_pre_iter,
                           cs_cdo_blas_dotprod_2pvsp,
                           cs_cdo_blas_square_norm_2pvsp);

  } /* Anderson acceleration */

  algo->prev_res = algo->res;

  double delta_pg = cs_cdo_blas_square_norm_pvsp_diff(pg_pre_iter, pg_cur_iter);
  double delta_pl = cs_cdo_blas_square_norm_pvsp_diff(pl_pre_iter, pl_cur_iter);

  algo->res = delta_pg + delta_pl;
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                possible divergence of the algorithm */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg_default(algo);

  if (algo->verbosity > 0) {

    if (algo->n_algo_iter == 1) {

      if (algo->verbosity > 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance"
                      "  ||D_Pg||  ||D_Pl||\n",
                      cs_param_get_nl_algo_label(nl_algo_type));
      else
        cs_log_printf(CS_LOG_DEFAULT,
                      "### GWF.TPF %10s.It    Algo.Res   Tolerance\n",
                      cs_param_get_nl_algo_label(nl_algo_type));

    }

    if (algo->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e %5.3e %5.3e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol,
                    sqrt(delta_pg), sqrt(delta_pl));
    else
      cs_log_printf(CS_LOG_DEFAULT,
                    "### GWF.TPF %10s.It%02d  %5.3e  %6.4e\n",
                    cs_param_get_nl_algo_label(nl_algo_type),
                    algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media
 *         with a non-linear resolution relying on the Picard algorithm and
 *         a coupled algorithm
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      model        type of hydraulic model (miscible/immiscible)
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_coupled_tpf_picard_compute(const cs_mesh_t              *mesh,
                            const cs_cdo_connect_t       *connect,
                            const cs_cdo_quantities_t    *cdoq,
                            const cs_time_step_t         *time_step,
                            cs_gwf_model_type_t           model,
                            cs_flag_t                     option_flag,
                            cs_gwf_tpf_t                 *mc)
{
  bool cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* A first resolution has been done followed by an update */

  cs_real_t  *pg_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pg_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  BFT_MALLOC(pg_k, 2*cdoq->n_vertices, cs_real_t);
  pl_k = pg_k + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val_pre, pg_k);
  cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val_pre, pl_k);

  pg_kp1 = mc->g_pressure->val;
  pl_kp1 = mc->l_pressure->val;

  /* Set the normalization factor */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_k);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_k);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  cs_log_printf(CS_LOG_DEFAULT, "%s: normalization=%6.4e\n",
                __func__, algo->normalization);

  /* Main non-linear loop */

  while (_check_nl_tpf_cvg(mc->nl_algo_type,
                           pg_k, pg_kp1, pl_k, pl_kp1,
                           algo) == CS_SLES_ITERATING) {

    cs_array_real_copy(cdoq->n_vertices, pg_kp1, pg_k);
    cs_array_real_copy(cdoq->n_vertices, pl_kp1, pl_k);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, mc->system);

    /* Get the soil state up to date with the last compute values for the
       liquid and gas pressures */

    /* relaxation */

    const double  relax = 0.;
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {

      pg_kp1[i] = (1 - relax)*pg_kp1[i] + relax*pg_k[i];
      pl_kp1[i] = (1 - relax)*pl_kp1[i] + relax*pl_k[i];

    }

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      model,
                      update_flag,
                      option_flag,
                      mc);

  } /* Until convergence */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          mc->system->param->name,
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pg_k);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media
 *         with a non-linear resolution relying on the Anderson algorithm and
 *         a coupled algorithm
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      model        type of hydraulic model (miscible/immiscible)
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_coupled_tpf_anderson_compute(const cs_mesh_t              *mesh,
                              const cs_cdo_connect_t       *connect,
                              const cs_cdo_quantities_t    *cdoq,
                              const cs_time_step_t         *time_step,
                              cs_gwf_model_type_t           model,
                              cs_flag_t                     option_flag,
                              cs_gwf_tpf_t                 *mc)
{
  bool  cur2prev = false;
  cs_flag_t  update_flag = 0;   /* No current to previous operation */
  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  /* A first resolution has been done followed by an update */

  cs_real_t  *pg_kp1 = NULL, *pl_kp1 = NULL;  /* at ^{n+1,k+1} */
  cs_real_t  *pg_k = NULL, *pl_k = NULL;      /* at ^{n+1,k} */

  BFT_MALLOC(pg_k, 2*cdoq->n_vertices, cs_real_t);
  pl_k = pg_k + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val_pre, pg_k);
  cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val_pre, pl_k);

  /* One needs only one array gathering the liquid and gas pressures */

  BFT_MALLOC(pg_kp1, 2*cdoq->n_vertices, cs_real_t);
  pl_kp1 = pg_kp1 + cdoq->n_vertices;

  cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val, pg_kp1);
  cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val, pl_kp1);

  /* Set the normalization factor */

  algo->normalization = cs_cdo_blas_square_norm_pvsp(pg_k);
  algo->normalization += cs_cdo_blas_square_norm_pvsp(pl_k);
  algo->normalization = sqrt(algo->normalization);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  /* Main non-linear loop */

  while (_check_nl_tpf_cvg(mc->nl_algo_type,
                           pg_k, pg_kp1, pl_k, pl_kp1,
                           algo) == CS_SLES_ITERATING) {

    cs_array_real_copy(cdoq->n_vertices, pg_kp1, pg_k);
    cs_array_real_copy(cdoq->n_vertices, pl_kp1, pl_k);

    /* Update the variables related to the groundwater flow system.
     * In case of an Anderson acceleration, pg_kp1 and pl_kp1 may be
     * updated */

    if (algo->n_algo_iter >= mc->anderson_param.starting_iter) {

      cs_array_real_copy(cdoq->n_vertices, pg_kp1, mc->g_pressure->val);
      cs_array_real_copy(cdoq->n_vertices, pl_kp1, mc->l_pressure->val);
    }

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      model,
                      update_flag,
                      option_flag,
                      mc);

    /* Build and solve the linear system related to the coupled system of
       equations. First call: current --> previous and then no operation */

    cs_equation_system_solve(cur2prev, mc->system);

    cs_array_real_copy(cdoq->n_vertices, mc->g_pressure->val, pg_kp1);
    cs_array_real_copy(cdoq->n_vertices, mc->l_pressure->val, pl_kp1);

  } /* Until convergence */

  /* Get the soil state up to date with the last compute values for the
     liquid and gas pressures */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    model,
                    update_flag,
                    option_flag,
                    mc);

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          mc->system->param->name,
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  /* Free temporary arrays and structures */

  BFT_FREE(pg_k);
  BFT_FREE(pg_kp1);
  cs_iter_algo_aa_free_arrays(algo->context);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new (unsteady) state for the groundwater flows module.
 *         Case of (miscible or immiscible) two-phase flows in porous media.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      model        type of hydraulic model (miscible/immiscible)
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

static void
_segregated_tpf_compute(const cs_mesh_t              *mesh,
                        const cs_time_step_t         *time_step,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *cdoq,
                        cs_gwf_model_type_t           model,
                        cs_flag_t                     option_flag,
                        cs_gwf_tpf_t                 *mc)
{
  assert(mc->use_incremental_solver);

  /* Copy the state for time t^n into field->val_pre */

  cs_field_current_to_previous(mc->g_pressure);
  cs_field_current_to_previous(mc->l_pressure);

  /* Since the cur2prev operation has been done, avoid to do it again */

  bool cur2prev = false;        /* Done just above */
  cs_flag_t  update_flag = 0;   /* No current to previous operation */

  /* Initialize the non-linear algorithm */

  cs_iter_algo_t  *algo = mc->nl_algo;
  assert(algo != NULL);

  cs_iter_algo_reset_nl(mc->nl_algo_type, algo);

  algo->normalization = cs_cdo_blas_square_norm_pvsp(mc->c_pressure->val);
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;
  else
    algo->normalization = sqrt(algo->normalization);

  /* First solve step:
   * 1. Solve the equation associated to Pl --> Pl^(n+1,1) knowing
   *    Pl^(n+1,0) = Pl^n, Pg^n and thus Pc^n and Sl^n
   *    --> Sl^(n+1,0) - Sl^n = 0 --> no source term
   * 2. Solve the equation associated to Pg --> Pg^(n+1,1) knowing
   *
   */

  cs_equation_solve(cur2prev, mesh, mc->wl_eq);

  cs_equation_solve(cur2prev, mesh, mc->hg_eq);

  /* Update the variables related to the groundwater flow system */

  cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                    model,
                    CS_FLAG_CURRENT_TO_PREVIOUS, /* Force this operation */
                    option_flag,
                    mc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
  char  label[32];

  sprintf(label, "Pl_iter%02d", algo->n_algo_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("liquid_pressure")->val,
                           time_step);

  sprintf(label, "Pg_iter%02d", algo->n_algo_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("gas_pressure")->val,
                           time_step);

  sprintf(label, "Pc_iter%02d", algo->n_algo_iter);
  cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                           CS_POST_WRITER_DEFAULT,
                           label,
                           1,
                           false,
                           false,
                           CS_POST_TYPE_cs_real_t,
                           cs_field_by_name("capillarity_pressure")->val,
                           time_step);
#endif

  cs_real_t  *dpc_kp1 = NULL;
  BFT_MALLOC(dpc_kp1, mesh->n_vertices, cs_real_t);

  _get_capillarity_pressure_increment(mesh, mc->wl_eq, mc->hg_eq, dpc_kp1);

  while(_check_nl_tpf_dpc_cvg(mc->nl_algo_type,
                              dpc_kp1, algo) == CS_SLES_ITERATING) {

    /* Solve step */

    cs_equation_solve(cur2prev, mesh, mc->wl_eq);

    cs_equation_solve(cur2prev, mesh, mc->hg_eq);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      model,
                      update_flag,
                      option_flag,
                      mc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_GWF_DBG > 1
    sprintf(label, "Pl_iter%02d", algo->n_algo_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("liquid_pressure")->val,
                             time_step);

    sprintf(label, "Pg_iter%02d", algo->n_algo_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("gas_pressure")->val,
                             time_step);

    sprintf(label, "Pc_iter%02d", algo->n_algo_iter);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_DEFAULT,
                             label,
                             1,
                             false,
                             false,
                             CS_POST_TYPE_cs_real_t,
                             cs_field_by_name("capillarity_pressure")->val,
                             time_step);
#endif

    _get_capillarity_pressure_increment(mesh, mc->wl_eq, mc->hg_eq, dpc_kp1);

  } /* while not converged */

  /* If something wrong happens, write a message in the listing */

  cs_iter_algo_post_check(__func__,
                          "Segregated incremental TPF solver",
                          cs_param_get_nl_algo_label(mc->nl_algo_type),
                          algo);

  BFT_FREE(dpc_kp1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize the model context structure for two-phase
 *        flows in a porous media
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tpf_t *
cs_gwf_tpf_create(void)
{
  cs_gwf_tpf_t  *mc = NULL;

  BFT_MALLOC(mc, 1, cs_gwf_tpf_t);

  /* Define the system of equations */
  /* ------------------------------ */

  /* Create a new equation associated to the mass conservation of water. The
     unknown is the capillarity pressure. This will stand for the (0,0)-block
     if a coupoled system is considered. */

  mc->wl_eq = cs_equation_add("w_conservation",        /* equation name */
                              "capillarity_pressure",  /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_PARAM_BC_HMG_NEUMANN);

  /* Create a new equation associated to the mass conservation of gas (hydrogen
     for instance). The unknown is the pressure in the gaseous phase. This will
     stand for the (1,1)-block if a coupled system is considered. */

  mc->hg_eq = cs_equation_add("h_conservation", /* equation name */
                              "gas_pressure",   /* variable name */
                              CS_EQUATION_TYPE_GROUNDWATER,
                              1,
                              CS_PARAM_BC_HMG_NEUMANN);

  /* Set to NULL the remaining pointers. Only a part is set according to the
     numerical settings */

  mc->wg_eqp = NULL;
  mc->hl_eqp = NULL;
  mc->system = NULL;

  /* Advection fields */
  /* ---------------- */

  /* Darcy flux (one assumes a CDOVB space discretization) */

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_GWF | CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  cs_adv_field_t  *l_adv_field = cs_advection_field_add("l_darcy_field",
                                                        adv_status);
  mc->l_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  mc->l_darcy->adv_field = l_adv_field;

  cs_adv_field_t  *g_adv_field = cs_advection_field_add("g_darcy_field",
                                                        adv_status);
  mc->g_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  mc->g_darcy->adv_field = g_adv_field;

  cs_adv_field_t  *t_adv_field = cs_advection_field_add("total_darcy_field",
                                                        adv_status);
  mc->t_darcy = cs_gwf_darcy_flux_create(cs_flag_dual_face_byc);
  mc->t_darcy->adv_field = t_adv_field;

  /* Properties
   * ----------
   *
   * Adding the properties related to the permeability (diffusion) is
   * postponed since one has to know which type of permeability is considered
   * (iso, ortho, or anisotropic) and this information is a result of the type
   * of soils which have been added.
   */

  /* Properties related to the water equation  */

  mc->time_wl_pty = NULL;
  mc->diff_wl_pty = NULL;
  mc->time_wg_pty = NULL;

  /* Properties related to the hydrogen equation  */

  mc->time_hg_pty = NULL;
  mc->diff_hg_pty = NULL;
  mc->reac_hg_pty = NULL;
  mc->time_hl_pty = NULL;
  mc->diff_hl_pty = NULL;

  /* Fields */
  /* ------ */

  mc->l_saturation = NULL;
  mc->c_pressure = NULL;
  mc->l_pressure = NULL;
  mc->g_pressure = NULL;

  /* Arrays */
  /* ------ */

  /* The properties will be defined using arrays.
   * Store these arrays of property values associated to equation terms */

  mc->time_wl_array = NULL;
  mc->diff_wl_array = NULL;
  mc->srct_wl_array = NULL;

  mc->time_wg_array = NULL;

  mc->time_hg_array = NULL;
  mc->diff_hg_array = NULL;
  mc->reac_hg_array = NULL;
  mc->srct_hg_array = NULL;

  mc->time_hl_array = NULL;
  mc->diff_hl_array = NULL;

  /* Array of additional variable or property values */

  mc->l_rel_permeability = NULL;
  mc->g_rel_permeability = NULL;
  mc->l_capacity = NULL;
  mc->capillarity_cell_pressure = NULL;
  mc->g_cell_pressure = NULL;
  mc->l_saturation_submesh = NULL;
  mc->l_saturation_submesh_pre = NULL;

  /* Model parameters (default values) */
  /* ---------------- */

  mc->l_mass_density = 1000;
  mc->l_viscosity = 1e-3;
  mc->g_viscosity = 2e-5;
  mc->l_diffusivity_h = 0;      /* immiscible case */
  mc->w_molar_mass = 18e-3;
  mc->h_molar_mass = 3e-3;
  mc->ref_temperature = 280;    /* in Kelvin */
  mc->henry_constant = 1e-20;   /* nearly immiscible case */

  /* Numerical parameters */

  mc->use_coupled_solver = false;
  mc->use_incremental_solver = true;
  mc->use_properties_on_submesh = true;
  mc->use_explicit_dsldt_liquid = false;

  mc->nl_algo_type = CS_PARAM_NL_ALGO_NONE; /* Linear algo. by default */

  mc->nl_algo_cvg.n_max_iter = 50;
  mc->nl_algo_cvg.rtol = 1e-5;
  mc->nl_algo_cvg.atol = 1e-10;
  mc->nl_algo_cvg.dtol = 1e3;

  mc->anderson_param.n_max_dir = 5;
  mc->anderson_param.starting_iter = 3;
  mc->anderson_param.max_cond = -1; /* No test by default */
  mc->anderson_param.beta = 1.0;    /* No damping by default */
  mc->anderson_param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

  mc->nl_algo = NULL;

  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the modelling of two-phase
 *        flows in a porous media
 *
 * \param[in, out] p_mc   pointer of pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_free(cs_gwf_tpf_t  **p_mc)
{
  if (p_mc == NULL)
    return;
  if (*p_mc == NULL)
    return;

  cs_gwf_tpf_t  *mc = *p_mc;

  /* System of equations are freed elsewhere (just after having freed
     cs_equation_t structures) */

  cs_gwf_darcy_flux_free(&(mc->l_darcy));
  cs_gwf_darcy_flux_free(&(mc->g_darcy));
  cs_gwf_darcy_flux_free(&(mc->t_darcy));

  BFT_FREE(mc->time_wl_array);
  BFT_FREE(mc->diff_wl_array);
  BFT_FREE(mc->srct_wl_array);
  BFT_FREE(mc->time_wg_array);
  BFT_FREE(mc->time_hg_array);
  BFT_FREE(mc->diff_hg_array);
  BFT_FREE(mc->reac_hg_array);
  BFT_FREE(mc->srct_hg_array);
  BFT_FREE(mc->time_hl_array);
  BFT_FREE(mc->diff_hl_array);

  BFT_FREE(mc->l_rel_permeability);
  BFT_FREE(mc->g_rel_permeability);
  BFT_FREE(mc->l_capacity);
  BFT_FREE(mc->capillarity_cell_pressure);
  BFT_FREE(mc->g_cell_pressure);
  BFT_FREE(mc->l_saturation_submesh);
  BFT_FREE(mc->l_saturation_submesh_pre);

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    /* Non-linearity: If the context is not NULL, this means that an Anderson
       algorithm has been activated otherwise nothing to do */

    cs_iter_algo_aa_free(mc->nl_algo);

    BFT_FREE(mc->nl_algo);

  }

  BFT_FREE(mc);
  *p_mc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup related to the model context of two-phase flows.
 *        Common to the different sub-models relying on two-phase flows.
 *
 * \param[in] model   model chosen for the hydraulic
 * \param[in] mc      pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_log_setup(cs_gwf_model_type_t    model,
                     cs_gwf_tpf_t          *mc)
{
  if (mc == NULL)
    return;

  if (model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
    _mtpf_log_setup(mc);
  else {
    assert(model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE);
    _itpf_log_setup(mc);
  }

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Reference temperature: %5.2f K\n",
                mc->ref_temperature);

  cs_gwf_darcy_flux_log(mc->l_darcy);
  cs_gwf_darcy_flux_log(mc->g_darcy);
  cs_gwf_darcy_flux_log(mc->t_darcy);

  if (mc->use_coupled_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Coupled solver\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Segregated solver\n");

  if (mc->use_incremental_solver)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Incremental solver\n");

  if (mc->use_properties_on_submesh)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Liquid saturation on submesh\n");

  if (mc->use_explicit_dsldt_liquid)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Explicit treatment dSl/dt (liq.)\n");

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Non-linear algo.: %s\n",
                  cs_param_get_nl_algo_name(mc->nl_algo_type));

    cs_log_printf(CS_LOG_SETUP, "  * GWF | Tolerances of non-linear algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e; max_iter: %d\n",
                  mc->nl_algo_cvg.rtol, mc->nl_algo_cvg.atol,
                  mc->nl_algo_cvg.dtol, mc->nl_algo_cvg.n_max_iter);

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      const cs_iter_algo_param_aa_t  aap = mc->anderson_param;

      cs_log_printf(CS_LOG_SETUP, "  * GWF | Anderson param: max. dir: %d; "
                    " start: %d; drop. tol: %5.3e; relax: %5.3e\n",
                    aap.n_max_dir, aap.starting_iter, aap.max_cond, aap.beta);
      cs_log_printf(CS_LOG_SETUP, "  * GWF | Anderson param: Dot product: %s\n",
                    cs_param_get_dotprod_type_name(aap.dp_type));

    }

  } /* There is a non-linear algorithm */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the model context according to the settings done inside
 *        the function cs_user_model()
 *        Case of a two-phase flows model in porous media
 *
 * \param[in, out] mc          pointer to the model context structure
 * \param[in, out] perm_type   type of permeability to handle
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init(cs_gwf_tpf_t            *mc,
                cs_property_type_t       perm_type)
{
  if (mc == NULL)
    return;

  if (mc->wl_eq == NULL || mc->hg_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  cs_equation_param_t  *wl_eqp = cs_equation_get_param(mc->wl_eq);
  cs_equation_param_t  *hg_eqp = cs_equation_get_param(mc->hg_eq);

  if (mc->use_coupled_solver) {

    /* Define the coupled system of equations */
    /* -------------------------------------- */

    /* Create the (0,1)-block related to the water in the gas phase */

    mc->wg_eqp = cs_equation_param_create("wg_block",
                                          CS_EQUATION_TYPE_GROUNDWATER,
                                          1,
                                          CS_PARAM_BC_HMG_NEUMANN);

    /* Create the (1,0)-block related to the hydrogen in the liquid phase */

    mc->hl_eqp = cs_equation_param_create("hl_block",
                                          CS_EQUATION_TYPE_GROUNDWATER,
                                          1,
                                          CS_PARAM_BC_HMG_NEUMANN);

    /* Add a 2x2 system of coupled equations and define each block */

    mc->system = cs_equation_system_add("TwoPhasePorousFlow",
                                        2,   /* system size */
                                        1);  /* scalar-valued block */

    /* Set all the blocks in the coupled system */

    cs_equation_system_assign_equation(0, mc->wl_eq, mc->system);  /* (0,0) */
    cs_equation_system_assign_equation(1, mc->hg_eq, mc->system);  /* (1,1) */
    cs_equation_system_assign_param(0, 1, mc->wg_eqp, mc->system); /* (0,1) */
    cs_equation_system_assign_param(1, 0, mc->hl_eqp, mc->system); /* (1,0) */

    /* Properties */
    /* ---------- */

    /* Properties which will be associated to blocks
     * - unsteady term for water eq. in the liquid phase      (0,0) block
     * - diffusion term for water eq. in the liquid phase     (0,0)-block
     * - diffusion term for water eq. in the gaseous phase    (0,1)-block
     * - unsteady term for hydrogen eq. in the gaseous phase  (1,1)-block
     * - unsteady term for hydrogen eq. in the liquid phase   (1,0)-block
     */

    mc->time_wl_pty = cs_property_add("time_wl_pty", CS_PROPERTY_ISO);
    mc->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
    mc->diff_wg_pty = cs_property_add("diff_wg_pty", perm_type);
    mc->time_hg_pty = cs_property_add("time_hg_pty", CS_PROPERTY_ISO);
    mc->time_hl_pty = cs_property_add("time_hl_pty", CS_PROPERTY_ISO);

    /* Add terms to the water equation */
    /* ------------------------------- */

    cs_equation_add_time(wl_eqp, mc->time_wl_pty);
    cs_equation_add_diffusion(wl_eqp, mc->diff_wl_pty);

    /* Cross-terms for the block (0,1) -- Water equation */

    cs_equation_add_diffusion(mc->wg_eqp, mc->diff_wg_pty);

    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_time(hg_eqp, mc->time_hg_pty);
    cs_equation_add_advection(hg_eqp, mc->t_darcy->adv_field);

    /* Cross-terms for the block (1,0) -- Hydrogen equation */

    cs_equation_add_time(mc->hl_eqp, mc->time_hl_pty);

  }
  else { /* Segregated solver */

    mc->use_incremental_solver = true; /* Segregated solver are always solved
                                          by increment */

    /* Properties */
    /* ---------- */

    mc->diff_wl_pty = cs_property_add("diff_wl_pty", perm_type);
    mc->diff_hg_pty = cs_property_add("diff_hg_pty", perm_type);

    if (mc->use_properties_on_submesh) {

      mc->time_hg_pty = cs_property_subcell_add("time_hg_pty", CS_PROPERTY_ISO);
      mc->reac_hg_pty = cs_property_subcell_add("reac_hg_pty", CS_PROPERTY_ISO);

    }
    else {

      mc->time_hg_pty = cs_property_add("time_hg_pty", CS_PROPERTY_ISO);
      mc->reac_hg_pty = cs_property_add("reac_hg_pty", CS_PROPERTY_ISO);

    }

    if (!mc->use_explicit_dsldt_liquid) {

      cs_property_t  *time_wl_pty = NULL;

      if (mc->use_properties_on_submesh)
        time_wl_pty = cs_property_subcell_add("time_wl_pty", CS_PROPERTY_ISO);
      else
        time_wl_pty = cs_property_add("time_wl_pty", CS_PROPERTY_ISO);

      mc->time_wl_pty = time_wl_pty;

    }

    /* Add terms to the water equation */
    /* ------------------------------- */

    if (!mc->use_explicit_dsldt_liquid)
      cs_equation_add_time(wl_eqp, mc->time_wl_pty);

    cs_equation_add_diffusion(wl_eqp, mc->diff_wl_pty);

    /* Add terms to the hydrogen equation */
    /* ---------------------------------- */

    cs_equation_add_time(hg_eqp, mc->time_hg_pty);
    cs_equation_add_diffusion(hg_eqp, mc->diff_hg_pty);
    cs_equation_add_reaction(hg_eqp, mc->reac_hg_pty);

  } /* Segregated or coupled solver */

  if (mc->use_incremental_solver) {

    wl_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;
    hg_eqp->incremental_algo_type = CS_PARAM_NL_ALGO_PICARD;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial setup stage for two-phase flows in porous media. At this
 *        stage, all soils have been defined and equation parameters are set.
 *        Case of a miscible or immiscible model.
 *
 * \param[in]      post_flag  optional postprocessing request(s)
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_setup(cs_flag_t         post_flag,
                      cs_gwf_tpf_t     *mc)
{
  if (mc == NULL)
    return;

  if (mc->wl_eq == NULL || mc->hg_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Equations are not defined for this model. Stop execution.\n",
              __func__);

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");
  const int  v_loc_id = cs_mesh_location_get_id_by_name("vertices");
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");

  /* Retrieve the pointers to the structure storing the equation parameters */

  cs_equation_param_t  *wl_eqp = cs_equation_get_param(mc->wl_eq);
  cs_equation_param_t  *hg_eqp = cs_equation_get_param(mc->hg_eq);

  assert(wl_eqp->space_scheme == hg_eqp->space_scheme);

  int loc_id = c_loc_id;
  if (wl_eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      wl_eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    loc_id = v_loc_id;
  else
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Add the variable fields (Keep always the previous state) */

  cs_equation_predefined_create_field(1, mc->wl_eq);
  cs_equation_predefined_create_field(1, mc->hg_eq);

  /* Set fields related to variables */

  mc->c_pressure = cs_equation_get_field(mc->wl_eq);
  mc->g_pressure = cs_equation_get_field(mc->hg_eq);

  /* One has to be consistent with the location of DoFs for the w_eq and h_eq
   * which are respectively related to the l_pressure and g_pressure */

  mc->l_pressure = cs_field_create("liquid_pressure",
                                   field_mask,
                                   loc_id,
                                   1,
                                   true); /* has_previous */

  cs_field_set_key_int(mc->l_pressure, log_key, 1);
  cs_field_set_key_int(mc->l_pressure, post_key, 1);

  /* Create a liquid saturation field attached to cells: S_l */

  int  pty_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

  mc->l_saturation = cs_field_create("liquid_saturation",
                                     pty_mask,
                                     c_loc_id,
                                     1,     /* dimension */
                                     true); /* has_previous */

  cs_field_set_key_int(mc->l_saturation, log_key, 1);
  if (post_flag & CS_GWF_POST_LIQUID_SATURATION)
    cs_field_set_key_int(mc->l_saturation, post_key, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]      flag      optional settings for the module
 * \param[in, out] mc        pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_finalize_setup(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq,
                          cs_flag_t                    flag,
                          cs_gwf_tpf_t                *mc)
{
  CS_NO_WARN_IF_UNUSED(flag);   /* will be useful for gravity effect */

  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  c2v_size = c2v->idx[n_cells];

  cs_equation_param_t  *wl_eqp = cs_equation_get_param(mc->wl_eq);
  cs_equation_param_t  *hg_eqp = cs_equation_get_param(mc->hg_eq);

  assert(cs_equation_get_type(mc->wl_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(cs_equation_get_type(mc->hg_eq) == CS_EQUATION_TYPE_GROUNDWATER);
  assert(wl_eqp->space_scheme == hg_eqp->space_scheme);

  if (wl_eqp->space_scheme != CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid space scheme", __func__);

  /* Set the Darcian flux (in the volume and at the boundary) */
  /* -------------------------------------------------------- */

  cs_gwf_darcy_flux_define(connect, cdoq, wl_eqp->space_scheme,
                           mc, _tpf_l_darcy_update,
                           mc->l_darcy);

  cs_gwf_darcy_flux_define(connect, cdoq, hg_eqp->space_scheme,
                           mc, _tpf_g_darcy_update,
                           mc->g_darcy);

  cs_gwf_darcy_flux_define(connect, cdoq, hg_eqp->space_scheme,
                           mc, _tpf_t_darcy_update,
                           mc->t_darcy);

  /* Allocate and initialize arrays for physical propoerties */
  /* ------------------------------------------------------- */

  /* Allocate and initialize the relative permeability in the liquid and gas
     phase */

  BFT_MALLOC(mc->l_rel_permeability, n_cells, cs_real_t);
  BFT_MALLOC(mc->g_rel_permeability, n_cells, cs_real_t);

  /* One assumes that the medium is saturated by default */

  cs_array_real_set_scalar(n_cells, 1., mc->l_rel_permeability);
  cs_array_real_set_scalar(n_cells, 1., mc->g_rel_permeability);

  /* Interpolation of the gas pressure at cell centers is used to define
   * the properties associated to terms in the water or hydrogen eq. */

  BFT_MALLOC(mc->g_cell_pressure, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->g_cell_pressure);

  if (mc->use_properties_on_submesh) {

    BFT_MALLOC(mc->l_capacity, c2v_size, cs_real_t);
    cs_array_real_fill_zero(c2v_size, mc->l_capacity);

  }
  else {

    BFT_MALLOC(mc->capillarity_cell_pressure, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->capillarity_cell_pressure);

    BFT_MALLOC(mc->l_capacity, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->l_capacity);

  }

  /* Properties related to terms in equations */
  /* ---------------------------------------- */

  /* Define the array storing the diffusion property for the water eq. */

  BFT_MALLOC(mc->diff_wl_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->diff_wl_array);

  cs_property_def_by_array(mc->diff_wl_pty,
                           NULL,                 /* all cells */
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_wl_array,
                           false,                /* not owner of the array */
                           true);                /* full length */

  /* Define the array storing the diffusion property in the hydrogen eq. */

  BFT_MALLOC(mc->diff_hg_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->diff_hg_array);

  cs_property_def_by_array(mc->diff_hg_pty,
                           NULL,                 /* all cells */
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hg_array,
                           false,                /* not owner of the array */
                           true);                /* full length */

  BFT_MALLOC(mc->diff_hl_array, n_cells, cs_real_t);
  cs_array_real_fill_zero(n_cells, mc->diff_hl_array);

  cs_property_def_by_array(mc->diff_hl_pty,
                           NULL,                 /* all cells */
                           cs_flag_primal_cell,  /* where data are located */
                           mc->diff_hl_array,
                           false,                /* not owner of the array */
                           true);                /* full length */

  if (mc->use_coupled_solver) {

    /* Define the arrays storing the time property for all blocks */

    BFT_MALLOC(mc->time_wl_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->time_wl_array);

    cs_property_def_by_array(mc->time_wl_pty,
                             NULL,                 /* all cells */
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_wl_array,
                             false,                /* not owner of the array */
                             true);                /* full length */

    BFT_MALLOC(mc->time_wg_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->time_wg_array);

    cs_property_def_by_array(mc->time_wg_pty,
                             NULL,                 /* all cells */
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_wg_array,
                             false,                /* not owner of the array */
                             true);                /* full length */

    BFT_MALLOC(mc->time_hl_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->time_hl_array);

    cs_property_def_by_array(mc->time_hl_pty,
                             NULL,
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_hl_array,
                             false,                /* not owner of the array */
                             true);                /* full length */

    BFT_MALLOC(mc->time_hg_array, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, mc->time_hg_array);

    cs_property_def_by_array(mc->time_hg_pty,
                             NULL,                 /* all cells */
                             cs_flag_primal_cell,  /* where data are located */
                             mc->time_hg_array,
                             false,                /* not owner of the array */
                             true);                /* full length */

  }
  else { /* Segregated system */

    assert(mc->use_incremental_solver);

    /* Treatment of the time for the water mass conservation
     * - Define the array storing the source term values (always)
     * - Define the array for the time property (only if implicit dSl/dt)
     */

    if (mc->use_properties_on_submesh) {

      BFT_MALLOC(mc->l_saturation_submesh, c2v_size, cs_real_t);
      cs_array_real_fill_zero(c2v_size, mc->l_saturation_submesh);

      BFT_MALLOC(mc->srct_wl_array, c2v_size, cs_real_t);
      cs_array_real_fill_zero(c2v_size, mc->srct_wl_array);

      cs_xdef_t  *st_def =
        cs_equation_add_source_term_by_array(cs_equation_get_param(mc->wl_eq),
                                             NULL,        /* all cells */
                                             cs_flag_dual_cell_byc,
                                             mc->srct_wl_array,
                                             false,       /* is owner ? */
                                             true);       /* full length */

      cs_xdef_array_set_adjacency(st_def, c2v);

      if (mc->use_explicit_dsldt_liquid) {

        BFT_MALLOC(mc->l_saturation_submesh_pre, c2v_size, cs_real_t);
        cs_array_real_fill_zero(c2v_size, mc->l_saturation_submesh_pre);

      }
      else { /* Implicit treatment for dSl/dt */

        BFT_MALLOC(mc->time_wl_array, c2v_size, cs_real_t);
        cs_array_real_fill_zero(c2v_size, mc->time_wl_array);

        cs_xdef_t  *pty_def =
          cs_property_def_by_array(mc->time_wl_pty,
                                   NULL,                   /* all cells */
                                   cs_flag_dual_cell_byc,  /* data location */
                                   mc->time_wl_array,
                                   false,                  /* not owner */
                                   true);                  /* full length */

        cs_xdef_array_set_adjacency(pty_def, c2v);

      }

    }
    else { /* Properties on cells */

      BFT_MALLOC(mc->srct_wl_array, n_cells, cs_real_t);
      cs_array_real_fill_zero(n_cells, mc->srct_wl_array);

      cs_equation_add_source_term_by_array(cs_equation_get_param(mc->wl_eq),
                                           NULL,   /* all cells */
                                           cs_flag_primal_cell,
                                           mc->srct_wl_array,
                                           false,  /* is owner ? */
                                           true);  /* full length */

      BFT_MALLOC(mc->time_wl_array, n_cells, cs_real_t);
      cs_array_real_fill_zero(n_cells, mc->time_wl_array);

      cs_property_def_by_array(mc->time_wl_pty,
                               NULL,                 /* all cells */
                               cs_flag_primal_cell,  /* data location */
                               mc->time_wl_array,
                               false,                /* not owner */
                               true);                /* full length */

    }

    /* Treatment of the source term for the hydrogen mass conservation */

    if (mc->use_properties_on_submesh) {

      BFT_MALLOC(mc->srct_hg_array, c2v_size, cs_real_t);
      cs_array_real_fill_zero(c2v_size, mc->srct_hg_array);

      cs_xdef_t  *st_def =
        cs_equation_add_source_term_by_array(cs_equation_get_param(mc->hg_eq),
                                             NULL,   /* all cells */
                                             cs_flag_dual_cell_byc,
                                             mc->srct_hg_array,
                                             false,  /* is owner ? */
                                             true);  /* full length */

      cs_xdef_array_set_adjacency(st_def, c2v);

    }
    else { /* Properties on cells */

      BFT_MALLOC(mc->srct_hg_array, n_cells, cs_real_t);
      cs_array_real_fill_zero(n_cells, mc->srct_hg_array);

      cs_equation_add_source_term_by_array(cs_equation_get_param(mc->hg_eq),
                                           NULL,   /* all cells */
                                           cs_flag_primal_cell,
                                           mc->srct_hg_array,
                                           false,  /* is owner ? */
                                           true);  /* full length */

    }

    /* Treatment of the time property for the hydrogen mass conservation */

    if (mc->use_properties_on_submesh) {

      BFT_MALLOC(mc->time_hg_array, c2v_size, cs_real_t);
      cs_array_real_fill_zero(c2v_size, mc->time_hg_array);

      cs_xdef_t  *pty_def =
        cs_property_def_by_array(mc->time_hg_pty,
                                 NULL,                  /* all cells */
                                 cs_flag_dual_cell_byc, /* data location */
                                 mc->time_hg_array,
                                 false,                 /* not owner */
                                 true);                 /* full length */

      cs_xdef_array_set_adjacency(pty_def, c2v);

    }
    else { /* Properties on cells */

      BFT_MALLOC(mc->time_hg_array, n_cells, cs_real_t);
      cs_array_real_fill_zero(n_cells, mc->time_hg_array);

      cs_property_def_by_array(mc->time_hg_pty,
                               NULL,                /* all cells */
                               cs_flag_primal_cell, /* where data are located */
                               mc->time_hg_array,
                               false,               /* not owner of the array */
                               true);               /* full length */

    }

    /* Treatment of the reaction property for the hydrogen mass conservation */

    if (mc->use_properties_on_submesh) {

      BFT_MALLOC(mc->reac_hg_array, c2v_size, cs_real_t);
      cs_array_real_fill_zero(c2v_size, mc->reac_hg_array);

      cs_xdef_t  *pty_def =
        cs_property_def_by_array(mc->reac_hg_pty,
                                 NULL,                  /* all cells */
                                 cs_flag_dual_cell_byc, /* data location */
                                 mc->reac_hg_array,
                                 false,                 /* not owner */
                                 true);                 /* full length */

      cs_xdef_array_set_adjacency(pty_def, c2v);

    }
    else {

      BFT_MALLOC(mc->reac_hg_array, n_cells, cs_real_t);
      cs_array_real_fill_zero(n_cells, mc->reac_hg_array);

      cs_property_def_by_array(mc->reac_hg_pty,
                               NULL,                /* all cells */
                               cs_flag_primal_cell, /* where data are located */
                               mc->reac_hg_array,
                               false,               /* not owner of the array */
                               true);               /* full length */
    }

  } /* Segregated solve */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last setup stage in the case of two-phase flows in a porous media
 *        (miscible or immiscible case)
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq        pointer to a cs_cdo_quantities_t structure
 * \param[in]      verbosity   level of berbosity for the gwf module
 * \param[in, out] mc          pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_init_values(const cs_cdo_connect_t        *connect,
                       const cs_cdo_quantities_t     *cdoq,
                       int                            verbosity,
                       cs_gwf_tpf_t                  *mc)
{
  if (mc == NULL)
    return;

  if (mc->use_incremental_solver) {

    /* This model is non-linear so that one needs an iterative algorithm to
       manage this non-linearity */

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_NONE)
      mc->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;

    /* Add a hook function to define the source term. This operation should
       be done after that the corresponding equation has been initialized
       (a builder structure has to be defined) */

    cs_equation_add_build_hook(mc->hg_eq,
                               mc,                 /* hook context */
                               _build_hg_incr_st); /* hook function */

    /* Initialize other quantities */

    if (mc->use_properties_on_submesh && mc->use_explicit_dsldt_liquid)
      cs_array_real_copy(connect->c2v->idx[connect->n_cells],
                         mc->l_saturation_submesh,       /* src */
                         mc->l_saturation_submesh_pre);  /* dest */

  }

  /* Initialize the non-linear algorithm if needed */

  if (mc->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    mc->nl_algo = cs_iter_algo_create(verbosity, mc->nl_algo_cvg);

    /* One assumes that the discretization schemes are all set to CDO
       vertex-based schemes */

    assert(mc->hg_eq->param->space_scheme == CS_SPACE_SCHEME_CDOVB);

    cs_lnum_t  size = cdoq->n_vertices;
    if (mc->use_coupled_solver)
      size *= 2;

    if (mc->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
      mc->nl_algo->context = cs_iter_algo_aa_create(mc->anderson_param,
                                                    size);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the new state for the groundwater flows module.
 *        Case of two-phase flows in porous media.
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a cs_time_step_t structure
 * \param[in]      model        type of hydraulic model (miscible/immiscible)
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the model context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_compute(const cs_mesh_t               *mesh,
                   const cs_cdo_connect_t        *connect,
                   const cs_cdo_quantities_t     *cdoq,
                   const cs_time_step_t          *time_step,
                   cs_gwf_model_type_t            model,
                   cs_flag_t                      option_flag,
                   cs_gwf_tpf_t                  *mc)
{
  if (mc == NULL)
    return;

  if (mc->use_coupled_solver) {

    bool cur2prev = true;

    /* Build and solve the linear system related to the coupled system of
     * equations. By default, a current to previous operation is performed so
     * that prev <-- n and cur <-- n+1,k=1 (the first guess for the next time
     * step) */

    cs_equation_system_solve(cur2prev, mc->system);

    /* Update the variables related to the groundwater flow system */

    cs_gwf_tpf_update(mesh, connect, cdoq, time_step,
                      model,
                      CS_FLAG_CURRENT_TO_PREVIOUS,
                      option_flag,
                      mc);

    switch (mc->nl_algo_type) {

    case CS_PARAM_NL_ALGO_PICARD:
      _coupled_tpf_picard_compute(mesh, connect, cdoq, time_step,
                                  model,
                                  option_flag,
                                  mc);
      break;

    case CS_PARAM_NL_ALGO_ANDERSON:
      _coupled_tpf_anderson_compute(mesh, connect, cdoq, time_step,
                                    model,
                                    option_flag,
                                    mc);
      break;

    default:
      break; /* Nothing else to do */
    }

  }
  else
    _segregated_tpf_compute(mesh, time_step, connect, cdoq,
                            model,
                            option_flag,
                            mc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the update step in the case of a two-phase flow model in
 *        porous media (miscible or immiscible)
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to a cs_cdo_quantities_t structure
 * \param[in]      ts           pointer to a cs_time_step_t structure
 * \param[in]      model        type of hydraulic model (miscible/immiscible)
 * \param[in]      update_flag  metadata associated to type of operation to do
 * \param[in]      option_flag  calculation option related to the GWF module
 * \param[in, out] mc           pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_tpf_update(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *cdoq,
                  const cs_time_step_t        *ts,
                  cs_gwf_model_type_t          model,
                  cs_flag_t                    update_flag,
                  cs_flag_t                    option_flag,
                  cs_gwf_tpf_t                *mc)
{
  CS_NO_WARN_IF_UNUSED(option_flag); /* For a later usage (gravity) */

  cs_real_t  time_eval = ts->t_cur;
  bool  cur2prev = false;

  if (update_flag & CS_FLAG_CURRENT_TO_PREVIOUS) {

    time_eval = cs_equation_get_time_eval(ts, mc->wl_eq);
    cur2prev = true;

  }

  if (cs_equation_get_space_scheme(mc->wl_eq) !=
      CS_SPACE_SCHEME_CDOVB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible space discretization.", __func__);

  /* Update the Darcy fluxes */

  cs_gwf_darcy_flux_t  *l_darcy = mc->l_darcy;

  l_darcy->update_func(connect,
                       cdoq,
                       mc->l_pressure->val,
                       time_eval,
                       cur2prev,
                       l_darcy);

  cs_gwf_darcy_flux_t  *g_darcy = mc->g_darcy;

  g_darcy->update_func(connect,
                       cdoq,
                       mc->g_pressure->val,
                       time_eval,
                       cur2prev,
                       g_darcy);

  cs_gwf_darcy_flux_t  *t_darcy = mc->t_darcy;

  t_darcy->update_func(connect,
                       cdoq,
                       NULL,    /* values are defined inside the call */
                       time_eval,
                       cur2prev,
                       t_darcy);

  /* New pressure values for the liquid and the gas have been computed (compute
     step) */

  const cs_real_t  *c_pr = mc->c_pressure->val;
  const cs_real_t  *g_pr = mc->g_pressure->val;

  /* Interpolate the vertex values to the cell values. The capillarity pressure
     at cells is the one used to update quantities related to a soil model */

  if (mc->capillarity_cell_pressure != NULL)
    cs_reco_pv_at_cell_centers(connect->c2v, cdoq, c_pr,
                               mc->capillarity_cell_pressure);

  /* Interpolate the pressure in the gaseous phase at cell centers */

  cs_reco_pv_at_cell_centers(connect->c2v, cdoq, g_pr,
                             mc->g_cell_pressure);

  /* Avoid to add an unsteady contribution at the first iteration  */

  if (update_flag & CS_FLAG_INITIALIZATION)
    cs_array_real_copy(connect->n_vertices, g_pr, mc->g_pressure->val_pre);

  /* Update the liquid pressure: P_l = P_g - P_c */

  if (cur2prev) {

    cs_field_current_to_previous(mc->l_pressure);
    cs_field_current_to_previous(mc->l_saturation);

    if (mc->use_properties_on_submesh && mc->use_explicit_dsldt_liquid)
      cs_array_real_copy(connect->c2v->idx[connect->n_cells],
                         mc->l_saturation_submesh,       /* src */
                         mc->l_saturation_submesh_pre);  /* dest */

  }

  cs_real_t  *l_pr = mc->l_pressure->val;

  /* Compute the value at vertices */

  for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++)
    l_pr[i] = g_pr[i] - c_pr[i];

  /* Update properties related to soils:
   * - l_rel_permeability, g_rel_permeability
   * - liquid_saturation
   * - capacity: \frac{\partial S_l}{\partial P_c}
   *
   * Either a call to a user-defined function or a predefined function if the
   * soil corresponds to a known model
   */

  cs_gwf_soil_update(time_eval, mesh, connect, cdoq);

  /* Define the liquid saturation in each cell when the liquid saturation has
     been defined on the submesh */

  if (mc->use_properties_on_submesh) {

    const cs_adjacency_t  *c2v = connect->c2v;
    cs_real_t  *l_sat = mc->l_saturation->val;

    /* Compute the average liquid saturation in each cell */

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      cs_real_t  _sat = 0;
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        _sat += mc->l_saturation_submesh[j] * cdoq->pvol_vc[j];
      l_sat[c_id] = _sat/cdoq->cell_vol[c_id];

    } /* Loop on cells */

    if (mc->use_explicit_dsldt_liquid &&
        (update_flag & CS_FLAG_INITIALIZATION))
      cs_array_real_copy(connect->c2v->idx[connect->n_cells],
                         mc->l_saturation_submesh,       /* src */
                         mc->l_saturation_submesh_pre);  /* dest */

  } /* Properties on submesh */

  /* Update arrays associated to property terms */

  int  dim = cs_property_get_dim(mc->diff_wl_pty);

  switch (dim) {

  case 1: /* Isotropic case */
    if (model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
      cs_gwf_soil_iso_update_mtpf_terms(mc);

    else { /* Immiscible model */

      if (mc->use_incremental_solver) {

        if (mc->use_properties_on_submesh)
          cs_gwf_soil_iso_update_itpf_terms_incr_submesh(ts, connect, mc);
        else
          cs_gwf_soil_iso_update_itpf_terms_incr(ts, mc);

      }
      else
        cs_gwf_soil_iso_update_itpf_terms(mc);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Only the isotropic case is available up to now.",
              __func__);
    break;

  } /* Switch on the permeability type */

  return time_eval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module in case
 *        of miscible or immiscible two-phase flows in porous media
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      post_flag  requested quantities to be postprocessed
 * \param[in, out] mc         pointer to the casted model context
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_op(const cs_cdo_connect_t          *connect,
                    const cs_cdo_quantities_t       *cdoq,
                    cs_flag_t                        post_flag,
                    cs_gwf_tpf_t                    *mc)
{
  assert(mc != NULL);

  if (cs_flag_test(post_flag, CS_GWF_POST_DARCY_FLUX_BALANCE) == false)
    return; /* Nothing to do */

  cs_gwf_darcy_flux_balance(connect, cdoq, cs_equation_get_param(mc->wl_eq),
                            mc->l_darcy);

  cs_gwf_darcy_flux_balance(connect, cdoq, cs_equation_get_param(mc->hg_eq),
                            mc->g_darcy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module
 *        in case of saturated two-phase flows (tpf) in porous media.
 *
 * \param[in] mesh_id      id of the output mesh for the current call
 * \param[in] n_cells      local number of cells of post_mesh
 * \param[in] cell_ids     list of cells (0 to n-1)
 * \param[in] post_flag    flag gathering quantities to postprocess
 * \param[in] abs_perm     property for the absolute permeability
 * \param[in] mc           pointer to the model context structure
 * \param[in] time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_tpf_extra_post(int                        mesh_id,
                      cs_lnum_t                  n_cells,
                      const cs_lnum_t            cell_ids[],
                      cs_flag_t                  post_flag,
                      const cs_property_t       *abs_perm,
                      const cs_gwf_tpf_t        *mc,
                      const cs_time_step_t      *time_step)
{
  if (mesh_id != CS_POST_MESH_VOLUME)
    return; /* Only postprocessings in the volume are defined */

  assert(mc != NULL);

  if (post_flag & CS_GWF_POST_SOIL_CAPACITY) {

    if (mc->l_capacity == NULL) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf("%s: Requested postprocessing for capacity but not"
                 " allocated.\n This postprocessing is skipped.\n", __func__);
    }
    else {

      if (mc->use_properties_on_submesh) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf("%s: Requested postprocessing for capacity at cells but"
                   " this functionnality is not yet available\n"
                   "%s: with CS_GWF_LIQUID_SATURATION_ON_SUBMESH\n",
                   __func__, __func__);
      }
      else
        cs_post_write_var(mesh_id,
                          CS_POST_WRITER_DEFAULT,
                          "l_capacity",
                          1,
                          false,  /* interlace */
                          true,   /* use_parent */
                          CS_POST_TYPE_cs_real_t,
                          mc->l_capacity,
                          NULL,
                          NULL,
                          time_step);

    }

  } /* Postprocess the soil capacity */

  if (post_flag & CS_GWF_POST_PERMEABILITY) {

    /* permeability = abs_permeability * l_rel_permeability */

    cs_real_t  *permeability = NULL;
    int  dim = cs_property_get_dim(abs_perm);
    int  post_dim = (dim == 1) ? 1 : 9;

    if (dim > 1) {

      BFT_MALLOC(permeability, post_dim*n_cells, cs_real_t);

      for (cs_lnum_t i = 0; i < n_cells; i++) {

        cs_real_t  tensor[3][3];

        cs_property_get_cell_tensor(cell_ids[i],
                                    time_step->t_cur,
                                    abs_perm,
                                    false, /* inversion */
                                    tensor);

        cs_real_t  *_cell_perm = permeability + post_dim*i;
        for (int ki = 0; ki < 3; ki++)
          for (int kj = 0; kj < 3; kj++)
            _cell_perm[3*ki+kj] = tensor[ki][kj];

      }

    }
    else {

      BFT_MALLOC(permeability, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        permeability[c_id] = cs_property_get_cell_value(cell_ids[c_id],
                                                        time_step->t_cur,
                                                        abs_perm);

    }

    assert(mc->l_rel_permeability != NULL);
    for (cs_lnum_t c = 0; c < n_cells; c++)
      for (int k = 0; k < post_dim; k++)
        permeability[post_dim*c+k] *= mc->l_rel_permeability[cell_ids[c]];

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "permeability",
                      post_dim,
                      true,     /* interlace */
                      false,    /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      permeability,
                      NULL,
                      NULL,
                      time_step);

    BFT_FREE(permeability);

  } /* Post-processing of the permeability field */

  if (post_flag & CS_GWF_POST_GAS_MASS_DENSITY) {

    const double  mh_ov_rt =
      mc->h_molar_mass / (mc->ref_temperature * cs_physical_constants_r);

    cs_real_t  *gas_mass_density = NULL;
    BFT_MALLOC(gas_mass_density, n_cells, cs_real_t);

    for (cs_lnum_t c = 0; c < n_cells; c++)
      gas_mass_density[c] = mh_ov_rt * mc->g_cell_pressure[c];

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_DEFAULT,
                      "gas_mass_density",
                      1,
                      true,
                      false,
                      CS_POST_TYPE_cs_real_t,
                      gas_mass_density,
                      NULL,
                      NULL,
                      time_step);

    BFT_FREE(gas_mass_density);

  } /* Post-processing of the gas mass density */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
