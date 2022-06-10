/*============================================================================
 * Handle the solidification module with CDO schemes
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdofb_scaleq.h"
#include "cs_navsto_system.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_solid_selection.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solidification.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solidification.c

  \brief Structure and functions handling the solidification module
         (modified Navier-Stokes + thermal module + transport equations)
*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_SOLIDIFICATION_DBG       0

static const char _state_names[CS_SOLIDIFICATION_N_STATES][32] = {

  "Solid",
  "Mushy",
  "Liquid",
  "Eutectic"

};

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static variables
 *============================================================================*/

static cs_solidification_t  *cs_solidification_structure = NULL;
static cs_real_t  cs_solidification_forcing_eps  = 1e-10;
static cs_real_t  cs_solidification_eutectic_threshold  = 1e-4;

static const double  cs_solidification_diffusion_eps = 1e-16;
static const char _err_empty_module[] =
  " Stop execution.\n"
  " The structure related to the solidification module is empty.\n"
  " Please check your settings.\n";

/*============================================================================
 * Private static inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the liquidus temperature knowing the bulk concentration
 *         Assumption of the lever rule.
 *
 * \param[in]  alloy    pointer to a binary alloy structure
 * \param[in]  conc     value of the bulk concentration
 *
 * \return the computed liquidus temperature
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_t_liquidus(const cs_solidification_binary_alloy_t     *alloy,
                const cs_real_t                             conc)
{
  return  fmax(alloy->t_eut, alloy->t_melt + alloy->ml * conc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the solidus temperature knowing the bulk concentration
 *         Assumption of the lever rule.
 *
 * \param[in]  alloy    pointer to a binary alloy structure
 * \param[in]  conc     value of the bulk concentration
 *
 * \return the computed solidus temperature
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_t_solidus(const cs_solidification_binary_alloy_t     *alloy,
               const cs_real_t                             conc)
{
  if (conc < alloy->cs1)
    return alloy->t_melt + alloy->ml * conc * alloy->inv_kp;
  else
    return alloy->t_eut;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of eta (Cliq = eta * Cbulk) knowing the bulk
 *         concentration and the phase diagram.
 *         Assumption of the lever rule.
 *
 * \param[in] alloy      pointer to a binary alloy structure
 * \param[in] conc       value of the bulk concentration
 *
 * \return the value of eta
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_eta(const cs_solidification_binary_alloy_t     *alloy,
         const cs_real_t                             conc)
{
  /* Update eta */

  if (conc > alloy->cs1)
    /* In this case Cl = C_eut = eta * Cbulk--> eta = C_eut/Cbulk */
    return alloy->c_eut/conc;
  else
    return alloy->inv_kp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine in which state is a couple (temp, conc)
 *         Assumption of the lever rule.
 *
 * \param[in]  alloy    pointer to a binary alloy structure
 * \param[in]  temp     value of the temperature
 * \param[in]  conc     value of the bulk concentration
 *
 * \return the state among (liquid, solid, mushy or eutectic)
 */
/*----------------------------------------------------------------------------*/

static inline cs_solidification_state_t
_which_state(const cs_solidification_binary_alloy_t     *alloy,
             const cs_real_t                             temp,
             const cs_real_t                             conc)
{
  const cs_real_t  t_liquidus = _get_t_liquidus(alloy, conc);

  if (temp > t_liquidus)
    return CS_SOLIDIFICATION_STATE_LIQUID;

  else {   /* temp < t_liquidus */

    const cs_real_t  t_solidus = _get_t_solidus(alloy, conc);
    if (temp > t_solidus)
      return CS_SOLIDIFICATION_STATE_MUSHY;

    else { /* temp < t_solidus */

      if (conc < alloy->cs1 || temp < alloy->t_eut_inf)
        return CS_SOLIDIFICATION_STATE_SOLID;
      else
        return CS_SOLIDIFICATION_STATE_EUTECTIC;

    } /* solidus */
  }   /* liquidus */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine in which state is a tuple (temp, conc, gl) from the
 *         evaluation of its enthalpy. The calling code has to be sure that the
 *         tuple is consistent.
 *         Assumption of the lever rule.
 *
 * \param[in]  alloy         pointer to a binary alloy structure
 * \param[in]  latent_heat   value of the latent heat coefficient
 * \param[in]  cp            value of the heat capacity
 * \param[in]  temp          value of the temperature
 * \param[in]  conc          value of the bulk concentration
 * \param[in]  gliq          value of the liquid fraction
 *
 * \return the state among (liquid, solid, mushy or eutectic)
 */
/*----------------------------------------------------------------------------*/

static inline cs_solidification_state_t
_which_state_by_enthalpy(const cs_solidification_binary_alloy_t    *alloy,
                         const cs_real_t                            latent_heat,
                         const cs_real_t                            cp,
                         const cs_real_t                            temp,
                         const cs_real_t                            conc,
                         const cs_real_t                            gliq)
{
  const cs_real_t  h_liq = cp*_get_t_liquidus(alloy, conc) + latent_heat;
  const cs_real_t  h = cp*temp + gliq*latent_heat;

  if (h > h_liq)
    return CS_SOLIDIFICATION_STATE_LIQUID;

  else {

    if (conc > alloy->cs1) {    /* Part with eutectic */

      const cs_real_t  h_sol = cp*alloy->t_eut;
      const cs_real_t  gleut = (conc - alloy->cs1)*alloy->dgldC_eut;
      const cs_real_t  h_eut = cp*alloy->t_eut + gleut*latent_heat;

      if (h > h_eut)
        return CS_SOLIDIFICATION_STATE_MUSHY;
      else if (h > h_sol)
        return CS_SOLIDIFICATION_STATE_EUTECTIC;
      else
        return CS_SOLIDIFICATION_STATE_SOLID;

    }
    else {                      /* Part without eutectic */

      const cs_real_t  h_sol = cp*(alloy->t_melt+alloy->ml*conc*alloy->inv_kp);
      if (h > h_sol)
        return CS_SOLIDIFICATION_STATE_MUSHY;
      else
        return CS_SOLIDIFICATION_STATE_SOLID;

    } /* Eutectic or not that is the question ? */

  } /* Liquid ? */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the derivatives of g_l w.r.t. the temperature and the
 *         bulk concentration when the current state is MUSHY
 *         Assumption of the lever rule.
 *
 * \param[in]  alloy    pointer to a binary alloy structure
 * \param[in]  temp     value of the temperature
 * \param[in]  conc     value of the bulk concentration
 * \param[out] dgldT    value of the derivative of g_l w.r.t. the temperature
 * \param[out] dgldC    value of the derivative of g_l w.r.t. the concentration
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_dgl_mushy(const cs_solidification_binary_alloy_t     *alloy,
               const cs_real_t                             temp,
               const cs_real_t                             conc,
               cs_real_t                                  *dgldT,
               cs_real_t                                  *dgldC)
{
  const double _dTm = temp - alloy->t_melt;
  const double _kml = alloy->ml * alloy->inv_kpm1;

  *dgldT =  _kml * conc/(_dTm*_dTm);
  *dgldC = -_kml / _dTm;
}

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the convergence of the non-linear algorithm
 *
 * \param[in]       nl_algo_type  type of non-linear algorithm
 * \param[in]       pre_iter      previous iterate values
 * \param[in, out]  cur_iter      current iterate values
 * \param[in, out]  algo          pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_check_nl_cvg(cs_param_nl_algo_t        nl_algo_type,
              const cs_real_t          *pre_iter,
              cs_real_t                *cur_iter,
              cs_iter_algo_t           *algo)
{
  assert(algo != NULL);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON && algo->n_algo_iter > 0) {

    /* TODO */

  } /* Anderson acceleration */

  algo->prev_res = algo->res;
  algo->res = cs_cdo_blas_square_norm_pcsp_diff(pre_iter, cur_iter);
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                divergence */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg(algo);

  if (algo->param.verbosity > 0) {

    if (algo->n_algo_iter == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### SOLIDIFICATION %12s.It      Algo.Res  Tolerance\n",
                    cs_param_get_nl_algo_label(nl_algo_type));
    cs_log_printf(CS_LOG_DEFAULT,
                  "### SOLIDIFICATION %12s.It%02d   %5.3e  %6.4e\n",
                  cs_param_get_nl_algo_label(nl_algo_type),
                  algo->n_algo_iter, algo->res, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the structure dedicated to the management of the
 *         solidification module
 *
 * \return a pointer to a new allocated cs_solidification_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_solidification_t *
_solidification_create(void)
{
  cs_solidification_t  *solid = NULL;

  BFT_MALLOC(solid, 1, cs_solidification_t);

  /* Default initialization */

  solid->model = CS_SOLIDIFICATION_N_MODELS;
  solid->options = 0;
  solid->post_flag = 0;
  solid->verbosity = 1;

  /* Properties */

  solid->mass_density = NULL;
  solid->viscosity = NULL;

  /* Quantities related to the liquid fraction */

  solid->g_l = NULL;
  solid->g_l_field = NULL;

  /* State related to each cell */

  solid->cell_state = NULL;

  /* Monitoring */

  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++)
    solid->n_g_cells[i] = 0;
  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++)
    solid->state_ratio[i] = 0;

  /* Plot writer related to the solidification process */

  solid->plot_state = NULL;

  /* Structure related to the thermal system solved as a sub-module */

  solid->temperature = NULL;
  solid->thermal_reaction_coef = NULL;
  solid->thermal_reaction_coef_array = NULL;
  solid->thermal_source_term_array = NULL;

  /* Structure cast on-the-fly w.r.t. the modelling choice */

  solid->model_context = NULL;

  /* Quantities/structure related to the forcing term treated as a reaction term
     in the momentum equation */

  solid->forcing_mom = NULL;
  solid->forcing_mom_array = NULL;
  solid->forcing_coef = 0;
  solid->first_cell = -1;

  return solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the number of cells in a given state for monitoring purpose
 *
 * \param[in]     connect    pointer to a cs_cdo_connect_t structure
 * \param[in]     quant      pointer to a cs_cdo_quantities_t structure
 * \param[in,out] solid      pointer to the main cs_solidification_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_monitor_cell_state(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant,
                    cs_solidification_t         *solid)
{
  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++) solid->n_g_cells[i] = 0;

  for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

    if (connect->cell_flag[c] & CS_FLAG_SOLID_CELL)
      solid->n_g_cells[CS_SOLIDIFICATION_STATE_SOLID] += 1;
    else
      solid->n_g_cells[solid->cell_state[c]] += 1;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the monitoring dedicated to the solidification module
 *
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_do_monitoring(const cs_cdo_quantities_t   *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  assert(solid->temperature != NULL);

  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++)
    solid->state_ratio[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  vol_c = quant->cell_vol[c_id];

    switch (solid->cell_state[c_id]) {
    case CS_SOLIDIFICATION_STATE_SOLID:
      solid->state_ratio[CS_SOLIDIFICATION_STATE_SOLID] += vol_c;
      break;
    case CS_SOLIDIFICATION_STATE_LIQUID:
      solid->state_ratio[CS_SOLIDIFICATION_STATE_LIQUID] += vol_c;
      break;
    case CS_SOLIDIFICATION_STATE_MUSHY:
      solid->state_ratio[CS_SOLIDIFICATION_STATE_MUSHY] += vol_c;
      break;
    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      solid->state_ratio[CS_SOLIDIFICATION_STATE_EUTECTIC] += vol_c;
      break;

    default: /* Should not be in this case */
      break;

    } /* End of switch */

  } /* Loop on cells */

  /* Finalize the monitoring step*/

  cs_parall_sum(CS_SOLIDIFICATION_N_STATES, CS_REAL_TYPE, solid->state_ratio);
  const double  inv_voltot = 100./quant->vol_tot;
  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++)
    solid->state_ratio[i] *= inv_voltot;

  cs_log_printf(CS_LOG_DEFAULT,
                "### Solidification monitoring: liquid/mushy/solid states\n"
                "  * Solid    | %6.2f\%% for %9lu cells;\n"
                "  * Mushy    | %6.2f\%% for %9lu cells;\n"
                "  * Liquid   | %6.2f\%% for %9lu cells;\n",
                solid->state_ratio[CS_SOLIDIFICATION_STATE_SOLID],
                solid->n_g_cells[CS_SOLIDIFICATION_STATE_SOLID],
                solid->state_ratio[CS_SOLIDIFICATION_STATE_MUSHY],
                solid->n_g_cells[CS_SOLIDIFICATION_STATE_MUSHY],
                solid->state_ratio[CS_SOLIDIFICATION_STATE_LIQUID],
                solid->n_g_cells[CS_SOLIDIFICATION_STATE_LIQUID]);

  if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
    cs_log_printf(CS_LOG_DEFAULT,
                  "  * Eutectic | %6.2f\%% for %9lu cells;\n",
                  solid->state_ratio[CS_SOLIDIFICATION_STATE_EUTECTIC],
                  solid->n_g_cells[CS_SOLIDIFICATION_STATE_EUTECTIC]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add a source term to the solute equation derived from an explicit
 *          use of the advective and diffusive operator
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
_fb_solute_source_term(const cs_equation_param_t     *eqp,
                       const cs_equation_builder_t   *eqb,
                       const void                    *eq_context,
                       const cs_cell_mesh_t          *cm,
                       void                          *context,
                       cs_hodge_t                    *mass_hodge,
                       cs_hodge_t                    *diff_hodge,
                       cs_cell_sys_t                 *csys,
                       cs_cell_builder_t             *cb)
{
  CS_UNUSED(context);
  CS_UNUSED(mass_hodge);
  CS_UNUSED(eqb);

  if (cb->cell_flag & CS_FLAG_SOLID_CELL)
    return; /* No solute evolution in permanent solid zone */

  const cs_cdofb_scaleq_t  *eqc = (const cs_cdofb_scaleq_t *)eq_context;

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  cs_real_t  *cl_c = alloy->c_l_cells;
  cs_real_t  *cl_f = alloy->c_l_faces;

  /* Diffusion part of the source term to add */

  cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                 diff_hodge);

  /* Define the local stiffness matrix: local matrix owned by the cellwise
     builder (store in cb->loc) */

  eqc->get_stiffness_matrix(cm, diff_hodge, cb);

  /* Build the cellwise array: c - c_l
     One should have c_l >= c. Therefore, one takes fmin(...,0) */

  for (short int f = 0; f < cm->n_fc; f++)
    cb->values[f] = fmin(csys->val_n[f] - cl_f[cm->f_ids[f]], 0);
  cb->values[cm->n_fc] = fmin(csys->val_n[cm->n_fc] - cl_c[cm->c_id], 0);

  /* Update the RHS with the diffusion contribution */

  cs_sdm_update_matvec(cb->loc, cb->values, csys->rhs);

  /* Define the local advection matrix */

  /* Open hook: Compute the advection flux for the numerical scheme and store
     the advection fluxes across primal faces */

  eqc->advection_open(eqp, cm, csys, eqc->advection_input, cb);

  eqc->advection_main(eqp, cm, csys, eqc->advection_scheme, cb);

  /* Build the cellwise array: c - c_l
     One should have c_l >= c. Therefore, one takes fmin(...,0) */

  for (short int f = 0; f < cm->n_fc; f++)
    cb->values[f] = fmin(csys->val_n[f] - cl_f[cm->f_ids[f]], 0);
  cb->values[cm->n_fc] = fmin(csys->val_n[cm->n_fc] - cl_c[cm->c_id], 0);

  /* Update the RHS with the convection contribution */

  cs_sdm_update_matvec(cb->loc, cb->values, csys->rhs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of (local) solid cells and enforce a zero-velocity
 *         for this selection
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_solid_cells(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  cs_gnum_t  n_solid_cells = solid->n_g_cells[CS_SOLIDIFICATION_STATE_SOLID];

  /* List of solid cells */

  cs_solid_selection_t  *scells = cs_solid_selection_get();

  if (n_solid_cells > (cs_gnum_t)scells->n_cells)
    BFT_REALLOC(scells->cell_ids, n_solid_cells, cs_lnum_t);

  scells->n_cells = n_solid_cells;

  if (n_solid_cells > 0) {

    n_solid_cells = 0;
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
      if (solid->cell_state[c_id] == CS_SOLIDIFICATION_STATE_SOLID)
        scells->cell_ids[n_solid_cells++] = c_id;
    }

    assert(n_solid_cells == solid->n_g_cells[CS_SOLIDIFICATION_STATE_SOLID]);

  }

  /* Parallel synchronization of the number of solid cells. Update the
     structure storing the list of solid cells and faces */

  cs_solid_selection_sync(connect);

  /* Enforce a zero velocity inside solid cells (enforcement of the momentum
     equation) */

  cs_navsto_system_set_solid_cells(scells->n_cells, scells->cell_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the enthalpy at each cell centers
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval      physical time at which evaluation is performed
 * \param[in]      temp        array of temperature values at each cell
 * \param[in]      g_l         array of the liquid fraction values at each cell
 * \param[in]      t_ref       reference temperature
 * \param[in]      latent_heat value of the latent heat coefficient
 * \param[in]      rho         property related to the mass density
 * \param[in]      cp          property related to the heat capacity
 * \param[in, out] enthalpy    array of enthalpy values at each cell
 */
/*----------------------------------------------------------------------------*/

static void
_compute_enthalpy(const cs_cdo_quantities_t    *quant,
                  cs_real_t                     t_eval,
                  const cs_real_t               temp[],
                  const cs_real_t               g_l[],
                  const cs_real_t               temp_ref,
                  const cs_real_t               latent_heat,
                  const cs_property_t          *rho,
                  const cs_property_t          *cp,
                  cs_real_t                     enthalpy[])
{
  assert(temp != NULL && g_l != NULL && enthalpy != NULL);

  if (quant->n_cells < 1)
    return;

  cs_real_t  rho_c, cp_c;

  bool  rho_is_uniform = cs_property_is_uniform(rho);
  bool  cp_is_uniform = cs_property_is_uniform(cp);

  /* Use cell with id 0 to evaluate the properties */

  if (rho_is_uniform)
    rho_c = cs_property_get_cell_value(0, t_eval, rho);

  if (cp_is_uniform)
    cp_c = cs_property_get_cell_value(0, t_eval, cp);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

    /* Retrieve the value of the properties if non uniform */

    if (!rho_is_uniform)
      rho_c = cs_property_get_cell_value(c, t_eval, rho);
    if (!cp_is_uniform)
      cp_c = cs_property_get_cell_value(c, t_eval, cp);

    enthalpy[c] = rho_c *
      /* part linked to the variation of  | part linked to the phase change
         temperature                      |                                 */
      ( cp_c * (temp[c] - temp_ref)       +  latent_heat * g_l[c] );

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*
 * Update functions for the Voller & Prakash modelling
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the liquid fraction, the cell state and the array
 *         used to compute the forcing term in the momentum equation.
 *         This corresponds to the methodology described in the paper
 *         written by Voller and Prakash (87).
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_voller_legacy(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;

  assert(solid->temperature != NULL);
  assert(v_model != NULL);

  cs_real_t  *g_l = solid->g_l_field->val;
  cs_real_t  *temp = solid->temperature->val;
  assert(temp != NULL);

  /* 1./(t_liquidus - t_solidus) = \partial g_l/\partial Temp */

  const cs_real_t  dgldT = 1./(v_model->t_liquidus - v_model->t_solidus);
  const cs_real_t  inv_forcing_eps = 1./cs_solidification_forcing_eps;

  assert(cs_property_is_uniform(solid->viscosity));
  const cs_real_t  viscl0 = cs_property_get_cell_value(solid->first_cell,
                                                       ts->t_cur,
                                                       solid->viscosity);
  const cs_real_t  forcing_coef = solid->forcing_coef * viscl0;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL) {

      g_l[c_id] = 0;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_SOLID;

    }

    /* Update the liquid fraction */

    else if (temp[c_id] < v_model->t_solidus) {

      g_l[c_id] = 0;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_SOLID;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */

      solid->forcing_mom_array[c_id] = forcing_coef*inv_forcing_eps;

    }
    else if (temp[c_id] > v_model->t_liquidus) {

      g_l[c_id] = 1;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_LIQUID;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */

      solid->forcing_mom_array[c_id] = 0;

    }
    else { /* Mushy zone */

      const cs_real_t  glc = (temp[c_id] - v_model->t_solidus) * dgldT;

      g_l[c_id] = glc;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_MUSHY;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */

      const cs_real_t  glm1 = 1 - glc;
      solid->forcing_mom_array[c_id] =
        forcing_coef * glm1*glm1/(glc*glc*glc + cs_solidification_forcing_eps);

    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the liquid fraction and the cell state. The array
 *         used to compute the forcing term in the momentum equation is not
 *         considered in the following function since there is no velocity
 *         field.  This corresponds to the methodology described in the paper
 *         written by Voller and Prakash (87).
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_voller_legacy_no_velocity(const cs_mesh_t               *mesh,
                                     const cs_cdo_connect_t        *connect,
                                     const cs_cdo_quantities_t     *quant,
                                     const cs_time_step_t          *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model =
    (cs_solidification_voller_t *)solid->model_context;

  assert(solid->temperature != NULL);
  assert(v_model != NULL);

  cs_real_t  *g_l = solid->g_l_field->val;
  cs_real_t  *temp = solid->temperature->val;
  assert(temp != NULL);

  /* 1./(t_liquidus - t_solidus) = \partial g_l/\partial Temp */

  const cs_real_t  dgldT = 1./(v_model->t_liquidus - v_model->t_solidus);

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL) {

      g_l[c_id] = 0;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_SOLID;

    }

    /* Update the liquid fraction */

    else if (temp[c_id] < v_model->t_solidus) {

      g_l[c_id] = 0;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_SOLID;

    }
    else if (temp[c_id] > v_model->t_liquidus) {

      g_l[c_id] = 1;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_LIQUID;

    }
    else { /* Mushy zone */

      const cs_real_t  glc = (temp[c_id] - v_model->t_solidus) * dgldT;

      g_l[c_id] = glc;
      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_MUSHY;

    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the reaction and source term for the thermal
 *         equation. This corresponds to the methodology described in the paper
 *         written by Voller and Prakash (87)
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_voller_legacy(const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;

  assert(v_model != NULL);
  assert(solid->temperature != NULL);

  const cs_real_t  *temp = solid->temperature->val;
  assert(temp != NULL);

  /* 1./(t_liquidus - t_solidus) = \partial g_l/\partial Temp */

  const cs_real_t  rho0 = solid->mass_density->ref_value;
  const cs_real_t  dgldT = 1./(v_model->t_liquidus - v_model->t_solidus);
  const cs_real_t  dgldT_coef = rho0*solid->latent_heat*dgldT/ts->dt[0];

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (solid->cell_state[c_id] == CS_SOLIDIFICATION_STATE_MUSHY) {

      solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
      solid->thermal_source_term_array[c_id] =
        dgldT_coef*temp[c_id]*quant->cell_vol[c_id];

    }
    else {

      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;

    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the reaction and source term for the thermal
 *         equation. One considers the state at the previous time step and that
 *         the kth sub-iteration to determine the solidification path and
 *         compute the related quantities.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_voller_path(const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;

  assert(v_model != NULL);
  assert(solid->temperature != NULL);

  const cs_real_t  *temp = solid->temperature->val;
  const cs_real_t  *temp_pre = solid->temperature->val_pre;
  assert(temp != NULL && temp_pre != NULL);

  /* 1./(t_liquidus - t_solidus) = \partial g_l/\partial Temp */

  const cs_real_t  dgldT = 1./(v_model->t_liquidus - v_model->t_solidus);
  const cs_real_t  coef =
    solid->mass_density->ref_value*solid->latent_heat/ts->dt[0];
  const cs_real_t  dgldT_coef = coef * dgldT;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL) {

      /* Keep the solid state during all the computation */

      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;

    }
    else if (temp[c_id] < v_model->t_solidus) {

      if (temp_pre[c_id] > v_model->t_liquidus) {

        /* Liquid --> solid state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*v_model->t_liquidus*quant->cell_vol[c_id];

      }
      else if (temp_pre[c_id] < v_model->t_solidus) {

        /* Solid --> Solid state */

        solid->thermal_reaction_coef_array[c_id] = 0;
        solid->thermal_source_term_array[c_id] = 0;

      }
      else { /* Mushy --> solid state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*temp_pre[c_id]*quant->cell_vol[c_id];

        /* Strictly speaking this should not be divided by 1/dt but with a
           smaller time step (Tsolidus is reached before the end of the time
           step) */

      }

    }
    else if (temp[c_id] > v_model->t_liquidus) {

      if (temp_pre[c_id] > v_model->t_liquidus) {

        /* Liquid --> liquid state */

        solid->thermal_reaction_coef_array[c_id] = 0;
        solid->thermal_source_term_array[c_id] = 0;

      }
      else if (temp_pre[c_id] < v_model->t_solidus) {

        /* Solid --> liquid state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*v_model->t_solidus*quant->cell_vol[c_id];

      }
      else { /* Mushy --> liquid state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*temp_pre[c_id]*quant->cell_vol[c_id];

      }

    }
    else {

      if (temp_pre[c_id] > v_model->t_liquidus) {

        /* Liquid --> mushy state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*v_model->t_liquidus*quant->cell_vol[c_id];

      }
      else if (temp_pre[c_id] < v_model->t_solidus) {

        /* Solid --> mushy state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*v_model->t_solidus*quant->cell_vol[c_id];

      }
      else { /* Mushy --> mushy state */

        solid->thermal_reaction_coef_array[c_id] = dgldT_coef;
        solid->thermal_source_term_array[c_id] =
          dgldT_coef*temp_pre[c_id]*quant->cell_vol[c_id];

      } /* State for the previous temp (n-1) */

    } /* State for the current temp (n+1,k+1) */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*
 * Update functions for the binary alloy modelling
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the state associated to each cell in the case of a binary
 *         alloy. No MPI synchronization has to be performed at this stage.
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_binary_alloy_final_state(const cs_cdo_connect_t      *connect,
                                 const cs_cdo_quantities_t   *quant,
                                 const cs_time_step_t        *ts)
{
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  /* Update the cell state (at this stage, one should have converged between
   * the couple (temp, conc) and the liquid fraction */

  const cs_real_t  *t_bulk = solid->temperature->val;
  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *g_l = solid->g_l_field->val;

  for (int i = 0; i < CS_SOLIDIFICATION_N_STATES; i++) solid->n_g_cells[i] = 0;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL) {

      solid->cell_state[c_id] = CS_SOLIDIFICATION_STATE_SOLID;
      solid->n_g_cells[CS_SOLIDIFICATION_STATE_SOLID] += 1;

    }
    else {

      cs_solidification_state_t
        state = _which_state_by_enthalpy(alloy,
                                         solid->latent_heat,
                                         solid->cp->ref_value,
                                         t_bulk[c_id],
                                         c_bulk[c_id],
                                         g_l[c_id]);

      solid->cell_state[c_id] = state;
      solid->n_g_cells[state] += 1;

    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the Darcy term (acting as a penalization) in the momentum
 *         equation and enforce solid cells by setting a zero mass flux.
 *         The parallel reduction on the cell state is performed here (not
 *         before to avoid calling the enforcement if no solid cell is locally
 *         detected).
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_velocity_forcing(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);

  cs_solidification_t  *solid = cs_solidification_structure;

  /* At this stage, the number of solid cells is a local count
   * Set the enforcement of the velocity for solid cells */

  _enforce_solid_cells(connect, quant);

  /* Parallel synchronization of the number of cells in each state
   * This should be done done now to avoid going to the cell enforcement whereas
   * there is nothing to do locally */

  cs_parall_sum(CS_SOLIDIFICATION_N_STATES, CS_GNUM_TYPE, solid->n_g_cells);

  assert(cs_property_is_uniform(solid->viscosity));
  const cs_real_t  viscl0 = cs_property_get_cell_value(0, ts->t_cur,
                                                       solid->viscosity);
  const cs_real_t  forcing_coef = solid->forcing_coef * viscl0;
  const cs_real_t  *g_l = solid->g_l_field->val;

  /* Set the forcing term in the momentum equation */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (g_l[c_id] < 1.) {       /* Not fully liquid */

      const cs_real_t gsc = 1 - g_l[c_id];
      const cs_real_t glc3 = g_l[c_id]*g_l[c_id]*g_l[c_id];

      solid->forcing_mom_array[c_id] =
        forcing_coef * gsc*gsc/(glc3 + cs_solidification_forcing_eps);

    }
    else
      solid->forcing_mom_array[c_id] = 0;

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the concentration of solute in the liquid phase at the cell
 *         center. This value is used in the buoyancy term in the momentum
 *         equation.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_clc(const cs_mesh_t             *mesh,
            const cs_cdo_connect_t      *connect,
            const cs_cdo_quantities_t   *quant,
            const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *t_bulk = solid->temperature->val;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;

  cs_real_t  *c_l = alloy->c_l_cells;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL) {
      c_l[c_id] = 0.;
      continue;
    }

    const cs_real_t  conc = c_bulk[c_id];
    const cs_real_t  temp = t_bulk[c_id];

    switch (_which_state(alloy, temp, conc)) {

    case CS_SOLIDIFICATION_STATE_SOLID:
      /* If this is the first time that one reaches the solid state for this
       * cell (i.e previously with g_l > 0), then one updates the liquid
       * concentration and one keeps that value */

      if (g_l_pre[c_id] > 0) {
        if (conc < alloy->cs1)
          c_l[c_id] = conc * alloy->inv_kp;
        else
          c_l[c_id] = alloy->c_eut;
      }
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      c_l[c_id] = (temp - alloy->t_melt) * alloy->inv_ml;
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      c_l[c_id] = conc;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      c_l[c_id] = alloy->c_eut;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the liquid fraction in each cell
 *         This function reproduces the same process as the one used in the
 *         legacy FV scheme.
 *         This corresponds to the case of a binary alloy model with no
 *         advective source term for the solute transport.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_legacy(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *t_bulk = solid->temperature->val;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;
  cs_real_t        *g_l = solid->g_l_field->val;

  /* Update g_l values in each cell as well as the cell state and the related
     count */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t  eta_new, gliq;

    const cs_real_t  eta_old = alloy->eta_coef_array[c_id];
    const cs_real_t  conc = c_bulk[c_id];
    const cs_real_t  temp = t_bulk[c_id];

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    /* Knowing in which part of the phase diagram we are and we then update the
     * value of the liquid fraction: g_l and eta (the coefficient between the
     * concentration in the liquid phase and the bulk concentration */

    switch (_which_state(alloy, temp, conc)) {

    case CS_SOLIDIFICATION_STATE_SOLID:
      gliq = 0.;
      if (g_l_pre[c_id] > 0)    /* Not in a solid state */
        eta_new = _get_eta(alloy, conc);
      else
        eta_new = eta_old;
     break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      gliq = alloy->inv_kpm1* (alloy->kp - alloy->ml*conc/(temp-alloy->t_melt));

      /* Make sure that the liquid fraction remains inside physical bounds */

      gliq = fmin(fmax(0, gliq), 1.);

      eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      gliq = 1;
      eta_new = 1;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      gliq = (conc - alloy->cs1)*alloy->dgldC_eut;

      /* Make sure that the liquid fraction remains inside physical bounds */

      gliq = fmin(fmax(0, gliq), 1.);

      eta_new = _get_eta(alloy, conc);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

    /* Update the liquid fraction and apply if needed a relaxation */

    if (alloy->gliq_relax > 0)
      g_l[c_id] = (1 - alloy->gliq_relax)*gliq + alloy->gliq_relax*g_l[c_id];
    else
      g_l[c_id] = gliq;

    /* Update eta and apply if needed a relaxation */

    if (alloy->eta_relax > 0)
      alloy->eta_coef_array[c_id] =
        (1-alloy->eta_relax)*eta_new + alloy->eta_relax*eta_old;
    else
      alloy->eta_coef_array[c_id] = eta_new;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the liquid fraction in each cell
 *         This function reproduces the same process as the one used in the
 *         legacy FV scheme.
 *         This corresponds to the case of a binary alloy model with an
 *         advective source term for the solute transport.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_legacy_ast(const cs_mesh_t             *mesh,
                      const cs_cdo_connect_t      *connect,
                      const cs_cdo_quantities_t   *quant,
                      const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *t_bulk = solid->temperature->val;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;
  cs_real_t        *g_l = solid->g_l_field->val;
  cs_real_t        *c_l = alloy->c_l_cells;

  /* Update g_l values in each cell as well as the cell state and the related
     count */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    cs_real_t  gliq = 1;        /* Initialization as liquid */

    const cs_real_t  conc = c_bulk[c_id];
    const cs_real_t  temp = t_bulk[c_id];

    /* Knowing in which part of the phase diagram we are and we then update
     * the value of the liquid fraction: g_l and the concentration of the
     * liquid "solute" c_l */

    switch (_which_state(alloy, temp, conc)) {

    case CS_SOLIDIFICATION_STATE_SOLID:
      gliq = 0.;

      /* If this is the first time that one reaches the solid state for this
       * cell (i.e previously with g_l > 0), then one updates the liquid
       * concentration and one keeps that value */

      if (g_l_pre[c_id] > 0) {
        if (conc < alloy->cs1)
          c_l[c_id] = conc * alloy->inv_kp;
        else
          c_l[c_id] = alloy->c_eut;
      }
     break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      gliq = alloy->inv_kpm1* (alloy->kp - alloy->ml*conc/(temp-alloy->t_melt));
      c_l[c_id] = (temp - alloy->t_melt) * alloy->inv_ml;
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      c_l[c_id] = conc;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      gliq = (conc - alloy->cs1)*alloy->dgldC_eut;
      c_l[c_id] = alloy->c_eut;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

    /* Make sure that the liquid fraction remains inside physical bounds */

    gliq = fmin(fmax(0, gliq), 1.);

    /* Relaxation if needed for the liquid fraction */

    if (alloy->gliq_relax > 0)
      g_l[c_id] = (1 - alloy->gliq_relax)*gliq + alloy->gliq_relax*g_l[c_id];
    else
      g_l[c_id] = gliq;

  } /* Loop on cells */

  /* Update c_l at face values */

  const cs_equation_t  *tr_eq = alloy->solute_equation;
  const cs_real_t  *c_bulk_f = cs_equation_get_face_values(tr_eq, false);
  const cs_real_t  *t_bulk_f = alloy->temp_faces;

  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++) {

    const cs_real_t  conc = c_bulk_f[f_id];
    const cs_real_t  temp = t_bulk_f[f_id];

    /* Knowing in which part of the phase diagram we are, we then update
     * the value of the concentration of the liquid "solute" */

    switch (_which_state(alloy, temp, conc)) {

    case CS_SOLIDIFICATION_STATE_SOLID:
      if (conc < alloy->cs1)
        alloy->c_l_faces[f_id] = conc * alloy->inv_kp;
      else
        alloy->c_l_faces[f_id] = alloy->c_eut;
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      alloy->c_l_faces[f_id] = (temp - alloy->t_melt) * alloy->inv_ml;
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      alloy->c_l_faces[f_id] = conc;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      alloy->c_l_faces[f_id] = alloy->c_eut;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for face %ld\n", __func__, (long)f_id);
      break;

    } /* Switch on face state */

  } /* Loop on faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the source term for the thermal equation.
 *         This function reproduces the same process as the one used in the
 *         legacy FV scheme. This corresponds to the binary alloy model.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_legacy(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *c_bulk_pre = alloy->c_bulk->val_pre;
  const cs_real_t  *t_bulk_pre = solid->temperature->val_pre;

  const cs_real_t  rhoL = solid->mass_density->ref_value * solid->latent_heat;
  const cs_real_t  rhoLovdt = rhoL/ts->dt[0];

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update: 0 by default */

    const cs_real_t  conc = c_bulk[c_id];
    const cs_real_t  conc_pre = c_bulk_pre[c_id];
    const cs_real_t  temp_pre = t_bulk_pre[c_id];

    /* Knowing in which part of the phase diagram we are, then we update
     * the value of the concentration of the liquid "solute" */

    switch (_which_state(alloy, temp_pre, conc_pre)) {

    case CS_SOLIDIFICATION_STATE_SOLID:
    case CS_SOLIDIFICATION_STATE_LIQUID:
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      {
        cs_real_t  dgldC, dgldT;
        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
        solid->thermal_source_term_array[c_id] =
          quant->cell_vol[c_id] * rhoLovdt * ( dgldT * temp_pre +
                                               dgldC * (conc_pre - conc) );
      }
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = quant->cell_vol[c_id] *
        rhoLovdt * alloy->dgldC_eut * (conc_pre - conc);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the liquid fraction in each cell and related quantities.
 *         This corresponds to the case of a binary alloy model with no
 *         advective source term for the solute transport.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_taylor(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const double  cpovL = solid->cp->ref_value/solid->latent_heat;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *c_bulk_pre = alloy->c_bulk->val_pre;
  const cs_real_t  *t_bulk_pre = solid->temperature->val_pre;
  cs_real_t        *t_bulk = solid->temperature->val;
  cs_real_t        *g_l = solid->g_l_field->val;

  /* Update g_l values in each cell as well as the cell state and the related
     count */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    const cs_real_t  conc = c_bulk[c_id];            /* conc_{n+1}^{k+1}  */
    const cs_real_t  temp = t_bulk[c_id];            /* temp_{n+1}^{k+1} */
    const cs_real_t  conc_pre = c_bulk_pre[c_id];
    const cs_real_t  temp_pre = t_bulk_pre[c_id];

    cs_real_t  dgldC, dgldT, gliq, eta_new;
    cs_solidification_state_t  state, state_pre;

    /* gliq, temp and conc iterates may not be related with the gliq(temp, conc)
     * function until convergence is reached. So one needs to be careful. */

    state = _which_state(alloy, temp, conc);
    state_pre = _which_state(alloy, temp_pre, conc_pre);
    eta_new = alloy->eta_coef_array[c_id]; /* avoid a warning */
    gliq = g_l[c_id];                      /* avoid a warning */

    /* Knowing in which part of the phase diagram we are and we then update
     * the value of the liquid fraction: g_l and the concentration of the
     * liquid "solute" */

    switch (state) {

    case CS_SOLIDIFICATION_STATE_SOLID:

      if (state_pre == CS_SOLIDIFICATION_STATE_LIQUID) {

        /* Liquid --> Solid transition */

        const cs_real_t  t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        const cs_real_t  t_star =
          ( cpovL*temp + dgldT*t_liquidus + dgldC*(conc_pre - conc) ) /
          ( cpovL + dgldT );

        t_bulk[c_id] = t_star;

        gliq = 1 + (dgldT*(t_star - t_liquidus) + dgldC*(conc-conc_pre));

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        if (t_star > alloy->t_eut_sup)  /* Mushy or liquid */
          eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
        else  /* Eutectic or solid */
          eta_new = _get_eta(alloy, conc);

      }
      else {
        gliq = 0;
        eta_new = _get_eta(alloy, conc);
      }
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      if (state_pre == CS_SOLIDIFICATION_STATE_LIQUID) {
        /* Liquid --> Mushy transition */
        const cs_real_t  t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        const cs_real_t  t_star =
          ( cpovL*temp + dgldT*t_liquidus + dgldC*(conc_pre - conc) ) /
          ( cpovL + dgldT );

        gliq = 1 + (dgldT*(t_star - t_liquidus) + dgldC*(conc-conc_pre));

        t_bulk[c_id] = t_star;

      }
      else
        gliq = alloy->inv_kpm1 *
          ( alloy->kp - alloy->ml*conc / (temp - alloy->t_melt) );

      /* Make sure that the liquid fraction remains inside physical bounds */

      gliq = fmin(fmax(0, gliq), 1.);

      eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      gliq = 1;
      eta_new = 1;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      if (state_pre == CS_SOLIDIFICATION_STATE_LIQUID) {

        /* Liquid --> Eutectic transition */

        const cs_real_t  t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        const cs_real_t  t_star =
          ( cpovL*temp + dgldT*t_liquidus + dgldC*(conc_pre - conc) ) /
          ( cpovL + dgldT );

        t_bulk[c_id] = t_star;

        gliq = 1 + (dgldT*(t_star - t_liquidus) + dgldC*(conc-conc_pre));

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        if (t_star > alloy->t_eut_inf)
          eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
        else
          eta_new = _get_eta(alloy, conc);

      }
      else {

        const cs_real_t  temp_k = alloy->tk_bulk[c_id];  /* temp_{n+1}^k */

        /* g_l[c_id] is the value at the iterate k */

        gliq = g_l[c_id] + cpovL * (temp_k - alloy->t_eut);

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        /* In this case Cl = C_eut = eta * Cbulk--> eta = C_eut/Cbulk */

        eta_new = _get_eta(alloy, conc);

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid state for cell %ld\n",
                __func__, (long)c_id);
      break;

    } /* Switch on cell state */

    /* Update the liquid fraction and apply if needed a relaxation */

    if (alloy->gliq_relax > 0)
      g_l[c_id] = (1 - alloy->gliq_relax)*gliq + alloy->gliq_relax*g_l[c_id];
    else
      g_l[c_id] = gliq;

    /* Update eta and apply if needed a relaxation */

    if (alloy->eta_relax > 0) {
      cs_real_t  eta_old = alloy->eta_coef_array[c_id];
      alloy->eta_coef_array[c_id] =
        (1-alloy->eta_relax)*eta_new + alloy->eta_relax*eta_old;
    }
    else
      alloy->eta_coef_array[c_id] = eta_new;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the source term for the thermal equation.
 *         This corresponds to the binary alloy model.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_taylor(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  cs_real_t  dgldC, dgldT;
  cs_solidification_state_t  state_k;

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *c_bulk_pre = alloy->c_bulk->val_pre;
  const cs_real_t  *t_bulk_pre = solid->temperature->val_pre;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;

  const cs_real_t  rhoL = solid->mass_density->ref_value * solid->latent_heat;
  const cs_real_t  rhoLovdt = rhoL/ts->dt[0];
  const double  cpovL = solid->cp->ref_value/solid->latent_heat;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    const cs_real_t  conc = c_bulk[c_id];
    const cs_real_t  conc_pre = c_bulk_pre[c_id];
    const cs_real_t  temp_pre = t_bulk_pre[c_id];
    const cs_real_t  gliq_pre = g_l_pre[c_id];

    const cs_real_t  rhocvolLovdt = quant->cell_vol[c_id] * rhoLovdt;

    state_k = _which_state(alloy, alloy->tk_bulk[c_id], alloy->ck_bulk[c_id]);

    /* Knowing in which part of the phase diagram we are, then we update
     * the value of the concentration of the liquid "solute" */

    switch (_which_state(alloy, temp_pre, conc_pre)) {

    case CS_SOLIDIFICATION_STATE_LIQUID:
      /* From the knowledge of the previous iteration, try something
         smarter... */

      if (state_k == CS_SOLIDIFICATION_STATE_LIQUID) {
        solid->thermal_reaction_coef_array[c_id] = 0;
        solid->thermal_source_term_array[c_id] = 0;
      }
      else { /* Liquid --> Mushy transition */
        /*      Liquid --> Solid transition */
        /*      Liquid --> Eutectic transition */

        const cs_real_t  t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
        solid->thermal_source_term_array[c_id] = rhocvolLovdt *
          ( dgldT * t_liquidus + dgldC * (conc_pre - conc) );

      }
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      {
        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
        solid->thermal_source_term_array[c_id] = rhocvolLovdt *
          ( dgldT * temp_pre + dgldC * (conc_pre - conc) );
      }
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      {
        const cs_real_t  temp_k = alloy->tk_bulk[c_id];  /* temp_{n+1}^k */

        solid->thermal_reaction_coef_array[c_id] = 0;

        /* Estimate the variation of liquid fraction so that the physical
           bounds are satisfied for the liquid fraction) */

        cs_real_t  dgl = cpovL * (temp_k - alloy->t_eut);

        if (dgl + gliq_pre < 0)
          dgl = -gliq_pre;
        else if (dgl + gliq_pre > 1)
          dgl = 1 - gliq_pre;

        solid->thermal_source_term_array[c_id] = rhocvolLovdt * dgl;

      }
      break;

    case CS_SOLIDIFICATION_STATE_SOLID:
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the liquid fraction in each cell and related quantities.
 *         This corresponds to the case of a binary alloy model with no
 *         advective source term for the solute transport.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_binary_path(const cs_mesh_t             *mesh,
                       const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  CS_UNUSED(ts);
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const double  L = solid->latent_heat;
  const cs_real_t  cp0 = solid->cp->ref_value;
  const double  cpovL = cp0/L;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *c_bulk_pre = alloy->c_bulk->val_pre;
  cs_real_t        *t_bulk = solid->temperature->val;
  const cs_real_t  *t_bulk_pre = solid->temperature->val_pre;
  cs_real_t        *g_l = solid->g_l_field->val;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;

  /* Update g_l values in each cell as well as the cell state and the related
     count */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    const cs_real_t  conc = c_bulk[c_id];            /* conc_{n+1}^{k+1}  */
    const cs_real_t  temp = t_bulk[c_id];            /* temp_{n+1}^{k+1} */
    const cs_real_t  conc_pre = c_bulk_pre[c_id];
    const cs_real_t  temp_pre = t_bulk_pre[c_id];
    const cs_real_t  gliq_pre = g_l_pre[c_id];

    cs_real_t  dgldC, dgldT, gliq, eta_new, t_liquidus, t_solidus;
    cs_real_t  c_star, t_star, dh, dgl;
    cs_solidification_state_t  state, state_pre;

    /* gliq, temp and conc iterates may not be related with the
     * gliq(temp, conc) function until convergence is reached. So one needs to
     * be careful. */

    gliq = gliq_pre; /* default initialization to avoid a warning */

    state = _which_state(alloy, temp, conc);
    state_pre = _which_state(alloy, temp_pre, conc_pre);
    eta_new = alloy->eta_coef_array[c_id];

    /* Knowing in which part of the phase diagram we are and we then update
     * the value of the liquid fraction: g_l and the concentration of the
     * liquid "solute" */

    switch (state) {

    case CS_SOLIDIFICATION_STATE_SOLID:
      /* ============================= */

      switch (state_pre) {
      case CS_SOLIDIFICATION_STATE_LIQUID: /* Liquid --> Solid transition */

        t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        t_star = ( cpovL*temp + 1 + dgldT*t_liquidus + dgldC*(conc_pre-conc) ) /
          ( cpovL + dgldT );

        gliq = 1 + (dgldT*(t_star - t_liquidus) + dgldC*(conc-conc_pre));

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        if (gliq > 0) {

          t_solidus = _get_t_solidus(alloy, conc);
          if (t_star > t_solidus) /* Mushy or liquid */
            eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );

          else {

            /* Remain on the solidus line and redefine a new state */

            t_star = t_solidus;
            eta_new = _get_eta(alloy, conc);
          }
        }
        else
          eta_new = _get_eta(alloy, conc);

        t_bulk[c_id] = t_star;
        break;

      case CS_SOLIDIFICATION_STATE_MUSHY: /* Mushy --> Solid transition */
        t_solidus = _get_t_solidus(alloy, conc);
        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        /* Variation of enthalpy when considering a mushy zone */

        dh = cp0*(temp - temp_pre) +
          L*(dgldC*(conc-conc_pre) + dgldT*(temp-temp_pre));

        if (conc < alloy->cs1) { /* without eutectic */

          /* Take into account the fact that the variation of gliq is only in
             the mushy zone */

          c_star = conc_pre +
            (dh - cp0*(temp-temp_pre) - dgldT*(t_solidus-temp_pre) )
            / (L*dgldC);

          gliq = gliq_pre + dgldT*(temp-temp_pre) + dgldC*(c_star-conc_pre);

          /* Make sure that the gliq remains inside physical bounds */

          gliq = fmin(fmax(0, gliq), 1.);
          if (gliq > 0) {        /* still in the mushy zone */
            eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
            t_bulk[c_id] = t_solidus + 1e-6;
          }
          else
            eta_new = _get_eta(alloy, conc);

        }
        else {                  /* with eutectic */

          c_star = conc +
            (dh - cp0*(t_solidus-temp_pre)
             - L*(dgldC*(conc-conc_pre) + dgldT*(t_solidus-temp_pre)) )
            / (L*alloy->dgldC_eut);

          if (c_star < alloy->cs1 || c_star > alloy->c_eut) {
            gliq = 0;
            eta_new = _get_eta(alloy, conc);
          }
          else {

            gliq = gliq_pre +
              dgldC*(conc-conc_pre) + dgldT*(t_solidus-temp_pre) +
              alloy->dgldC_eut * (c_star - conc);

            /* Make sure that the gliq remains inside physical bounds */

            gliq = fmin(fmax(0, gliq), 1.);
            if (gliq > 0)         /* remains on the eutectic plateau */
              t_bulk[c_id] = t_solidus;

            eta_new = _get_eta(alloy, c_star);

          } /* Invalid c_star */

        } /* Eutectic transition taken into account */
        break;

      case CS_SOLIDIFICATION_STATE_EUTECTIC: /* Eutectic --> Solid transition */
        _get_dgl_mushy(alloy, alloy->t_eut, conc_pre, &dgldT, &dgldC);

        /* Variation of gl when considering how is implemented the eutectic
           zone */

        dgl = dgldT*(temp-temp_pre) + alloy->dgldC_eut*(conc-conc_pre);
        dh = cp0*(temp -temp_pre) + dgl*L;

        /* If one remains on the eutectic plateau, then the concentration should
           be c_star w.r.t. dh = dgldC_eut * (C* - Cn) since Tk+1 = Tn = Teut */

        c_star = conc_pre + dh/(L*alloy->dgldC_eut);

        if (c_star < alloy->cs1 || c_star > alloy->c_eut) {

          /* In fact the final state is below the eutectic plateau */

          gliq = 0;
          eta_new = _get_eta(alloy, conc);

        }
        else {

          gliq = gliq_pre + alloy->dgldC_eut*(c_star-conc_pre);
          gliq = fmin(fmax(0, gliq), 1.);
          eta_new = _get_eta(alloy, c_star);
          if (gliq > 0)
            t_bulk[c_id] = alloy->t_eut;

        }
        break;

      default: /* Solid --> solid */
        gliq = 0;
        if (gliq_pre > 0) /* Otherwise keep the same value for eta */
          eta_new = _get_eta(alloy, conc);
        break;

      } /* Switch on the previous state */
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      /* ============================= */

      switch (state_pre) {

      case CS_SOLIDIFICATION_STATE_LIQUID: /* Liquid --> Mushy transition */
        t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        t_star = ( cpovL*temp + dgldT*t_liquidus + dgldC*(conc_pre - conc) ) /
          ( cpovL + dgldT );

        gliq = 1 + (dgldT*(t_star-t_liquidus) + dgldC*(conc-conc_pre));

        t_bulk[c_id] = t_star;
        break;

      case CS_SOLIDIFICATION_STATE_MUSHY:
        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        gliq = gliq_pre +
          (dgldT*(temp-temp_pre) + dgldC*(conc-conc_pre));
        break;

      default:
        gliq = alloy->inv_kpm1 *
          ( alloy->kp - alloy->ml*conc / (temp - alloy->t_melt) );

      } /* End of switch on the previous state */

      /* Make sure that the liquid fraction remains inside physical bounds */
      gliq = fmin(fmax(0, gliq), 1.);

      eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
      break;

    case CS_SOLIDIFICATION_STATE_LIQUID:
      /* ============================== */

      gliq = 1;
      eta_new = 1;
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      /* ================================ */

      switch (state_pre) {

      case CS_SOLIDIFICATION_STATE_LIQUID: /* Liquid --> Eutectic transition */
        t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        t_star = ( cpovL*temp + dgldT*t_liquidus + dgldC*(conc_pre - conc) ) /
          ( cpovL + dgldT );

        t_bulk[c_id] = t_star;

        gliq = 1 + (dgldT*(t_star - t_liquidus) + dgldC*(conc-conc_pre));

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        if (t_star > alloy->t_eut_inf)
          eta_new = 1/( gliq * (1-alloy->kp) + alloy->kp );
        else
          eta_new = _get_eta(alloy, conc);
        break;

      case CS_SOLIDIFICATION_STATE_MUSHY: /* Mushy --> Eutectic transition */
        assert(conc > alloy->cs1);

        /* First part of the path in the mushy zone */

        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        gliq = g_l_pre[c_id] +
          alloy->dgldC_eut*(conc - conc_pre) + dgldT*(alloy->t_eut - temp_pre);

        /* Make sure that the liquid fraction remains inside physical bounds */
        gliq = fmin(fmax(0, gliq), 1.);

        eta_new = _get_eta(alloy, conc);
        break;

      default: /* eutectic --> eutectic or solid --> eutectic */

        _get_dgl_mushy(alloy, alloy->t_eut, conc_pre, &dgldT, &dgldC);

        /* Variation of gl when considering how is implemented the eutectic
           zone */

        dgl = dgldT*(temp-temp_pre) + alloy->dgldC_eut*(conc-conc_pre);
        dh = cp0*(temp -temp_pre) + dgl*L;

        /* If one remains on the eutectic plateau, then the concentration
           should be c_star w.r.t. dh = dgldC_eut * (C* - Cn) since
           Tk+1 = Tn = Teut */

        c_star = conc_pre + dh/(L*alloy->dgldC_eut);

        if (c_star < alloy->cs1 || c_star > alloy->c_eut) {

          gliq = (conc - alloy->cs1)*alloy->dgldC_eut;

          /* In this case Cl = C_eut = eta * Cbulk--> eta = C_eut/Cbulk */

          eta_new = _get_eta(alloy, conc);

        }
        else {

          gliq = gliq_pre + alloy->dgldC_eut*(c_star-conc_pre);
          if (gliq > 0)         /* Remains on the eutectic plateau */
            t_bulk[c_id] = alloy->t_eut;

          eta_new = _get_eta(alloy, c_star);

        }

        /* Make sure that the liquid fraction remains inside physical bounds */

        gliq = fmin(fmax(0, gliq), 1.);

        break;

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid state for cell %ld\n",
                __func__, (long)c_id);
      break;

    } /* Switch on cell state */

    /* Update the liquid fraction and apply if needed a relaxation */

    if (alloy->gliq_relax > 0)
      g_l[c_id] = (1 - alloy->gliq_relax)*gliq + alloy->gliq_relax*g_l[c_id];
    else
      g_l[c_id] = gliq;

    /* Update eta and apply if needed a relaxation */

    if (alloy->eta_relax > 0) {
      cs_real_t  eta_old = alloy->eta_coef_array[c_id];
      alloy->eta_coef_array[c_id] =
        (1-alloy->eta_relax)*eta_new + alloy->eta_relax*eta_old;
    }
    else
      alloy->eta_coef_array[c_id] = eta_new;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the source term for the thermal equation.
 *         This corresponds to the binary alloy model.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_binary_path(const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *ts)
{
  CS_UNUSED(mesh);
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  const cs_real_t  *c_bulk = alloy->c_bulk->val;
  const cs_real_t  *c_bulk_pre = alloy->c_bulk->val_pre;
  const cs_real_t  *t_bulk = solid->temperature->val;
  const cs_real_t  *t_bulk_pre = solid->temperature->val_pre;

  const cs_real_t  rhoL = solid->mass_density->ref_value * solid->latent_heat;
  const cs_real_t  rhoLovdt = rhoL/ts->dt[0];

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
      continue; /* No update */

    const cs_real_t  conc_kp1 = c_bulk[c_id]; /* Solute transport solved */
    const cs_real_t  conc_k = alloy->ck_bulk[c_id];
    const cs_real_t  temp_k = t_bulk[c_id];

    const cs_real_t  conc_pre = c_bulk_pre[c_id];
    const cs_real_t  temp_pre = t_bulk_pre[c_id];

    const cs_real_t  rhocvolLovdt = quant->cell_vol[c_id] * rhoLovdt;
    cs_real_t  dgldC, dgldT, t_solidus, t_liquidus;

    cs_solidification_state_t  state_k = _which_state(alloy, temp_k, conc_k);

    /* Knowing in which part of the phase diagram we are, then we update the
     * value of the concentration of the liquid "solute" */

    switch (_which_state(alloy, temp_pre, conc_pre)) {

    case CS_SOLIDIFICATION_STATE_LIQUID:
      /* ==============================
       * From the knowledge of the previous iteration, try something smarter...
       */

      switch (state_k) {

      case CS_SOLIDIFICATION_STATE_MUSHY:    /* Liquid --> Mushy */
        t_liquidus = _get_t_liquidus(alloy, conc_pre);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
        solid->thermal_source_term_array[c_id] = rhocvolLovdt *
          ( dgldT * t_liquidus + dgldC * (conc_pre-conc_kp1) );
        break;

      case CS_SOLIDIFICATION_STATE_EUTECTIC: /* Liquid --> Eutectic */
      case CS_SOLIDIFICATION_STATE_SOLID:    /* Liquid --> Solid */
        t_liquidus = _get_t_liquidus(alloy, conc_pre);
        t_solidus = _get_t_solidus(alloy, conc_kp1);

        _get_dgl_mushy(alloy, t_liquidus, conc_pre, &dgldT, &dgldC);

        solid->thermal_reaction_coef_array[c_id] = 0;
        solid->thermal_source_term_array[c_id] = rhocvolLovdt *
          ( dgldT * (t_liquidus-t_solidus) + dgldC * (conc_pre-conc_kp1) );
        break;

      default:                 /* Liquid */
        solid->thermal_reaction_coef_array[c_id] = 0;
        solid->thermal_source_term_array[c_id] = 0;
        break;

      } /* End of switch on the state k */
      break;

    case CS_SOLIDIFICATION_STATE_MUSHY:
      /* ============================= */
      {
        _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);

        switch (state_k) {

        case CS_SOLIDIFICATION_STATE_SOLID:    /* Mushy --> Solid transition */
          solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
          if (conc_kp1 < alloy->cs1)  /* Part without eutectic */
            solid->thermal_source_term_array[c_id] = rhocvolLovdt *
              ( dgldT*temp_pre + dgldC*(conc_pre-conc_kp1) );
          else                        /* Part with eutectic */
            solid->thermal_source_term_array[c_id] = rhocvolLovdt *
              ( dgldT*temp_pre + alloy->dgldC_eut*(conc_pre-conc_kp1) );
          break;

        case CS_SOLIDIFICATION_STATE_EUTECTIC: /* Mushy --> Eutectic */
          assert(conc_kp1 > alloy->cs1);

          /* First part in the mushy zone but not the final part */

          solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
          solid->thermal_source_term_array[c_id] = rhocvolLovdt *
            ( dgldT*temp_pre + alloy->dgldC_eut*(conc_pre-conc_kp1) );
          break;

        default:
          solid->thermal_reaction_coef_array[c_id] = dgldT * rhoLovdt;
          solid->thermal_source_term_array[c_id] = rhocvolLovdt *
            ( dgldT * temp_pre + dgldC * (conc_pre - conc_kp1) );
          break;

        }

      }
      break;

    case CS_SOLIDIFICATION_STATE_EUTECTIC:
      /* ================================ */
      {
        cs_real_t  r_coef = 0;
        cs_real_t  s_coef = alloy->dgldC_eut * (conc_pre - conc_kp1);

        if (solid->options & CS_SOLIDIFICATION_WITH_PENALIZED_EUTECTIC) {
          if (state_k == CS_SOLIDIFICATION_STATE_EUTECTIC ||
              state_k == CS_SOLIDIFICATION_STATE_SOLID) {
            if (conc_kp1 > alloy->cs1 && conc_kp1 < alloy->c_eut) {
              _get_dgl_mushy(alloy, temp_pre, conc_pre, &dgldT, &dgldC);
              r_coef = dgldT * rhoLovdt;
              s_coef += dgldT*alloy->t_eut;
            }
          }
        }

        solid->thermal_reaction_coef_array[c_id] = r_coef;
        solid->thermal_source_term_array[c_id] = rhocvolLovdt * s_coef;

      }
      break;

    case CS_SOLIDIFICATION_STATE_SOLID:
      /* ============================= */
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid state for cell %ld\n", __func__, (long)c_id);
      break;

    } /* Switch on cell state */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the source term for the thermal equation.
 *         This function is related to the Stefan problem with a liquid fraction
 *         being a step function w.r.t. the temperature
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_thm_stefan(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts)
{
  if (mesh->n_cells < 1)
    return;

  cs_real_t  rho_c, rhoLovdt;

  cs_solidification_t  *solid = cs_solidification_structure;

  const cs_real_t  Lovdt = solid->latent_heat/ts->dt[0];
  const cs_real_t  *g_l = solid->g_l_field->val;
  const cs_real_t  *g_l_pre = solid->g_l_field->val_pre;
  const cs_real_t  *vol = quant->cell_vol;

  bool  rho_is_uniform = cs_property_is_uniform(solid->mass_density);

  /* Use the first cell to set the value */

  if (rho_is_uniform) {
    rho_c = cs_property_get_cell_value(0, ts->t_cur, solid->mass_density);
    rhoLovdt = rho_c * Lovdt;
  }

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

    /* Retrieve the value of the properties */

    if (!rho_is_uniform) {
      rho_c = cs_property_get_cell_value(c, ts->t_cur, solid->mass_density);
      rhoLovdt = rho_c * Lovdt;
    }

    if (connect->cell_flag[c] & CS_FLAG_SOLID_CELL) /* Tag as solid for all the
                                                       computation */
      continue; /* No update: 0 by default */

    /* reaction_coef_array is set to zero. Only the source term is updated */

    if (fabs(g_l[c] - g_l_pre[c]) > 0)
      solid->thermal_source_term_array[c] = rhoLovdt*vol[c]*(g_l_pre[c]-g_l[c]);
    else
      solid->thermal_source_term_array[c] = 0;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the liquid fraction in each cell and the temperature if
 *         needed. Case of the Stefan model.
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_gl_stefan(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  const cs_time_step_t        *ts)
{
  CS_UNUSED(ts);

  if (mesh->n_cells < 1)
    return;

  cs_real_t  cp_c, cpovL;

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_stefan_t  *model =
    (cs_solidification_stefan_t *)solid->model_context;

  bool  cp_is_uniform = cs_property_is_uniform(solid->cp);

  cs_real_t  *temp = solid->temperature->val;
  cs_real_t  *g_l = solid->g_l_field->val;

  /* Use the first cell to set the value */

  if (cp_is_uniform) {
    cp_c = cs_property_get_cell_value(0, ts->t_cur, solid->cp);
    cpovL = cp_c/solid->latent_heat;
  }

  /* Update g_l values in each cell as well as the cell state and the related
     count */

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

    /* Retrieve the value of the property */

    if (!cp_is_uniform) {
      cp_c = cs_property_get_cell_value(c, ts->t_cur, solid->cp);
      cpovL = cp_c/solid->latent_heat;
    }

    if (connect->cell_flag[c] & CS_FLAG_SOLID_CELL)
      continue;  /* Tag as solid during all the computation
                    => No update: 0 by default */

    if (temp[c] > model->t_change) {

      if (g_l[c] < 1) {  /* Not in a stable state */

        /* Compute a new g_l */

        g_l[c] += cpovL * (temp[c] - model->t_change);
        if (g_l[c] < 1) {

          temp[c] = model->t_change;
          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_MUSHY;

        }
        else { /* Overshoot of the liquid fraction */

          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_LIQUID;
          temp[c] = model->t_change + 1./cpovL * (g_l[c] - 1);
          g_l[c] = 1.;

        }

      }
      else {  /* g_l = 1, stable state */

        solid->cell_state[c] = CS_SOLIDIFICATION_STATE_LIQUID;

      }

    }
    else { /* T < T_ch */

      if (g_l[c] > 0) { /* Not in a stable state */

        /* Compute a new g_l */

        g_l[c] += cpovL * (temp[c] - model->t_change);

        if (g_l[c] < 0) {       /* Undershoot of the liquid fraction */

          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_SOLID;
          temp[c] = model->t_change + 1./cpovL * g_l[c];
          g_l[c] = 0.;

        }
        else {

          temp[c] = model->t_change;
          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_MUSHY;

        }

      }
      else { /* g_l = 0, stable state */

        solid->cell_state[c] = CS_SOLIDIFICATION_STATE_SOLID;

      }

    }

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function aims at computing the new couple
 *         (temperature,liquid fraction) defining the new state
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_stefan_thermal_non_linearities(const cs_mesh_t              *mesh,
                                const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *quant,
                                const cs_time_step_t         *time_step)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_stefan_t  *s_model
    = (cs_solidification_stefan_t *)solid->model_context;

  const size_t  csize = quant->n_cells*sizeof(cs_real_t);
  const cs_equation_t  *t_eq = solid->thermal_sys->thermal_eq;

  /* Retrieve the current values */

  cs_real_t  *temp = cs_equation_get_cell_values(t_eq, false);
  cs_real_t  *g_l = solid->g_l_field->val;
  cs_real_t  *enthalpy = solid->enthalpy->val;

  cs_real_t  *hk = NULL;        /* enthalpy h^{n+1,k} */
  BFT_MALLOC(hk, quant->n_cells, cs_real_t);
  memcpy(hk, enthalpy, csize);

  /* Non-linear iterations (k) are performed to converge on the relation
   * h^{n+1,k+1} = h^{n+1,k} + eps with eps a user-defined tolerance
   *
   * h = rho.cp (Temp - Tref) + rho.L.gliq
   *
   * One sets:
   * T^{n+1,0} = T^n and gl^{n+1,0} = gl^n
   */

  cs_equation_current_to_previous(t_eq);
  cs_field_current_to_previous(solid->g_l_field);
  cs_field_current_to_previous(solid->enthalpy);

  /* Initialize loop stopping criteria */

  cs_real_t  delta_h = 1 + s_model->max_delta_h;
  int iter = 0;

  while ( delta_h > s_model->max_delta_h && iter < s_model->n_iter_max) {

    /* Compute the new thermal source term */

    s_model->update_thm_st(mesh, connect, quant, time_step);

    /* Solve the thermal system */

    cs_thermal_system_compute(false, /* No cur2prev inside a non-linear
                                        iterative process */
                              mesh, connect, quant, time_step);

    /* Compute the new liquid fraction (and update the temperature if needed) */

    s_model->update_gl(mesh, connect, quant, time_step);

    /* Now compute the enthalpy knowing the temperature and the liquid
     * fraction.
     * enthalpy stores k+1,n+1 and hk stores k,n+1
     */

    _compute_enthalpy(quant,
                      time_step->t_cur,     /* t_eval */
                      temp,                 /* temperature */
                      g_l,                  /* liquid fraction */
                      s_model->t_change,    /* temp_ref */
                      solid->latent_heat,   /* latent heat coeff. */
                      solid->mass_density,  /* rho */
                      solid->cp,            /* cp */
                      enthalpy);            /* computed enthalpy */

    delta_h = -1;
    for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

      cs_real_t  dh = fabs(enthalpy[c] - hk[c]);
      hk[c] = enthalpy[c];

      if (dh > delta_h)
        delta_h = dh;

    } /* Loop on cells */

    iter++;
    if (solid->verbosity > 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "### Solidification.NL: k= %d | delta_enthalpy= %5.3e\n",
                    iter, delta_h);

  } /* Until convergence */

  BFT_FREE(hk);

  /* Monitoring */

  if (solid->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  "## Solidification: Stop after %d iters, delta = %5.3e\n",
                  iter, delta_h);

  _monitor_cell_state(connect, quant, solid);

  /* Parallel synchronization of the number of cells in each state */

  cs_parall_sum(CS_SOLIDIFICATION_N_STATES, CS_GNUM_TYPE, solid->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new system state (temperature, liquid fraction) using
 *         the methodology defined in Voller & Prakash (87)
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_voller_prakash_87(const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant,
                   const cs_time_step_t         *time_step)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Solidification process with a pure component without segregation */

  cs_solidification_voller_t  *v_model =
    (cs_solidification_voller_t *)solid->model_context;

  /* Solve the thermal system */

  cs_thermal_system_compute(true, /* operate a cur2prev operation inside */
                            mesh, connect, quant, time_step);

  /* Update fields and properties which are related to solved variables */

  cs_field_current_to_previous(solid->g_l_field);

  v_model->update_gl(mesh, connect, quant, time_step);

  v_model->update_thm_st(mesh, connect, quant, time_step);

  /* Post-processing */

  if (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY)
    _compute_enthalpy(quant,
                      time_step->t_cur,        /* t_eval */
                      solid->temperature->val, /* temperature */
                      solid->g_l_field->val,   /* liquid fraction */
                      v_model->t_solidus,      /* temp_ref */
                      solid->latent_heat,      /* latent heat coeff. */
                      solid->mass_density,     /* rho */
                      solid->cp,               /* cp */
                      solid->enthalpy->val);   /* computed enthalpy */

  /* Monitoring */

  _monitor_cell_state(connect, quant, solid);

  /* At this stage, the number of solid cells is a local count
   * Set the enforcement of the velocity for solid cells */

  _enforce_solid_cells(connect, quant);

  /* Parallel synchronization of the number of cells in each state (It should
     be done after _enforce_solid_cells() */

  cs_parall_sum(CS_SOLIDIFICATION_N_STATES, CS_GNUM_TYPE, solid->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the new system state (temperature, liquid fraction) using
 *         the methodology defined in Voller & Prakash (87) but also taking
 *         into account the non-linearities stemming from the thermal source
 *         term
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_voller_non_linearities(const cs_mesh_t              *mesh,
                        const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *quant,
                        const cs_time_step_t         *time_step)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Solidification process with a pure component without segregation */

  cs_solidification_voller_t  *v_model =
    (cs_solidification_voller_t *)solid->model_context;

  cs_iter_algo_t  *algo = v_model->nl_algo;

  assert(algo != NULL);

  const size_t  csize = quant->n_cells*sizeof(cs_real_t);

  /* Retrieve the current values */

  cs_real_t  *hkp1 = solid->enthalpy->val;
  cs_real_t  *hk = NULL;   /* enthalpy  ^{n+1,k} */
  BFT_MALLOC(hk, quant->n_cells, cs_real_t);

  /* Initialize the stopping criteria */

  cs_iter_algo_reset_nl(v_model->nl_algo_type, algo);

  algo->normalization = sqrt(cs_cdo_blas_square_norm_pcsp(hkp1));
  if (algo->normalization < cs_math_zero_threshold)
    algo->normalization = 1.0;

  /* Non-linear iterations (k) are performed to converge on the relation
   * h^{n+1,k+1} = h^{n+1,k} + eps with eps a user-defined tolerance
   *
   * h = rho.cp (Temp - Tref) + rho.L.gliq
   *
   * One sets:
   * T^{n+1,0} = T^n and gl^{n+1,0} = gl^n
   */

  cs_equation_current_to_previous(solid->thermal_sys->thermal_eq);
  cs_field_current_to_previous(solid->g_l_field);
  cs_field_current_to_previous(solid->enthalpy);

  do {

    /* Compute the new thermal source term */

    v_model->update_thm_st(mesh, connect, quant, time_step);

    /* Solve the thermal system */

    cs_thermal_system_compute(false, /* No cur2prev inside a non-linear
                                        iterative process */
                              mesh, connect, quant, time_step);

    /* Compute the new liquid fraction (and update the temperature if needed) */

    v_model->update_gl(mesh, connect, quant, time_step);

    /* Now compute the enthalpy knowing the temperature and the liquid
     * fraction.
     * enthalpy stores k+1,n+1 and hk stores k,n+1
     */

    memcpy(hk, hkp1, csize);

    _compute_enthalpy(quant,
                      time_step->t_cur,        /* t_eval */
                      solid->temperature->val, /* temperature */
                      solid->g_l_field->val,   /* liquid fraction */
                      v_model->t_solidus,      /* temp_ref */
                      solid->latent_heat,      /* latent heat coeff. */
                      solid->mass_density,     /* rho */
                      solid->cp,               /* cp */
                      hkp1);                   /* computed enthalpy */

  } /* Until convergence */
  while (_check_nl_cvg(v_model->nl_algo_type,
                       hk, hkp1, algo) == CS_SLES_ITERATING);

  BFT_FREE(hk);

  /* Monitoring */

  if (solid->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "## Solidification: Stop non-linear algo. after %d iters,"
                  " residual = %5.3e\n",
                  algo->n_algo_iter, algo->res);

  /* Monitoring */

  _monitor_cell_state(connect, quant, solid);

  /* At this stage, the number of solid cells is a local count
   * Set the enforcement of the velocity for solid cells */

  _enforce_solid_cells(connect, quant);

  /* Parallel synchronization of the number of cells in each state (It should
     be done after _enforce_solid_cells() */

  cs_parall_sum(CS_SOLIDIFICATION_N_STATES, CS_GNUM_TYPE, solid->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function aims at computing the new temperature/bulk concentration
 *         state for the next iteration as well as updating all related
 *         quantities
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_default_binary_coupling(const cs_mesh_t              *mesh,
                         const cs_cdo_connect_t       *connect,
                         const cs_cdo_quantities_t    *quant,
                         const cs_time_step_t         *time_step)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;
  cs_equation_t  *c_eq = alloy->solute_equation;

  const size_t  csize = quant->n_cells*sizeof(cs_real_t);
  const cs_equation_t  *t_eq = solid->thermal_sys->thermal_eq;
  const cs_real_t  rho0 = solid->mass_density->ref_value;

  cs_real_t  *temp = cs_equation_get_cell_values(t_eq, false);
  cs_real_t  *conc = cs_equation_get_cell_values(c_eq, false);
  cs_real_t  *g_l = solid->g_l_field->val;

  /* Compute the state at t^(n+1) knowing that at state t^(n) */

  if (solid->options & CS_SOLIDIFICATION_USE_EXTRAPOLATION) {

    /* At this stage (i.e. before previous to current: val = n, val_pre = n-1 */

    cs_real_t  *temp_pre = cs_equation_get_cell_values(t_eq, true);
    cs_real_t  *conc_pre = cs_equation_get_cell_values(c_eq, true);

    /* Extrapolation at f_{n+1} = 2*f_n - f_{n-1} */

    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
      alloy->tx_bulk[c_id] = 2*temp[c_id] - temp_pre[c_id];
      alloy->cx_bulk[c_id] = 2*conc[c_id] - conc_pre[c_id];
    }

  } /* Extrapolation is requested */

  /* Non-linear iterations (k) are also performed to converge on the relation
   * gliq^{k+1} = gliq(temp^{k+1}, conc^{k+1})
   *
   * Cbulk^{0}_{n+1} = Cbulk_{n}
   * Tbulk^{0}_{n+1} = Tbulk_{n}
   * gl^{0}_{n+1} = gl_{n}
   */

  cs_equation_current_to_previous(c_eq);
  cs_equation_current_to_previous(t_eq);
  cs_field_current_to_previous(solid->g_l_field);

  /* At the beginning, field_{n+1}^{k=0} = field_n */

  memcpy(alloy->tk_bulk, temp, csize);
  memcpy(alloy->ck_bulk, conc, csize);

  cs_real_t  delta_temp = 1 + alloy->delta_tolerance;
  cs_real_t  delta_cbulk = 1 + alloy->delta_tolerance;

  alloy->iter = 0;
  while ( ( delta_temp  > alloy->delta_tolerance ||
            delta_cbulk > alloy->delta_tolerance  ) &&
          alloy->iter   < alloy->n_iter_max) {

    /* Solve Cbulk^(k+1)_{n+1} knowing Cbulk^{k}_{n+1}  */

    cs_equation_solve(false,  /* No cur2prev inside a non-linear iterative
                                 process */
                      mesh, c_eq);

    /* Update the source term for the thermal equation */

    alloy->update_thm_st(mesh, connect, quant, time_step);

    /* Solve the thermal system */

    cs_thermal_system_compute(false, /* No cur2prev inside a non-linear
                                        iterative process */
                              mesh, connect, quant, time_step);

    /* Update fields and properties which are related to solved variables
     * g_l, state */

    alloy->update_gl(mesh, connect, quant, time_step);

    /* Update the diffusion property related to the solute */

    if (alloy->diff_coef > cs_solidification_diffusion_eps) {

      const double  rho_D = rho0 * alloy->diff_coef;

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < quant->n_cells; i++)
        alloy->diff_pty_array[i] = (g_l[i] > 0) ?
          rho_D * g_l[i] : cs_solidification_diffusion_eps;

    }

    /* Evolution of the temperature and the bulk concentration during this
       iteration */

    delta_temp = -1, delta_cbulk = -1;
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      cs_real_t  dtemp = fabs(temp[c_id] - alloy->tk_bulk[c_id]);
      cs_real_t  dconc = fabs(conc[c_id] - alloy->ck_bulk[c_id]);

      alloy->tk_bulk[c_id] = temp[c_id];
      alloy->ck_bulk[c_id] = conc[c_id];

      if (dtemp > delta_temp)
        delta_temp = dtemp;
      if (dconc > delta_cbulk)
        delta_cbulk = dconc;

    } /* Loop on cells */

    /* Parallel synchronization */

    if (cs_glob_n_ranks > 1) {

      cs_real_t  parall_delta_max[2] = {delta_temp, delta_cbulk};

      cs_parall_max(2, CS_REAL_TYPE, parall_delta_max);

      delta_temp = parall_delta_max[0];
      delta_cbulk = parall_delta_max[1];

    }

    alloy->iter += 1;
    if (solid->verbosity > 0) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "### Solidification.NL: "
                    " k= %d | delta_temp= %5.3e | delta_cbulk= %5.3e\n",
                    alloy->iter, delta_temp, delta_cbulk);
    }

  } /* while iterating */

  /* Update the liquid concentration of the solute (c_l) */

  alloy->update_clc(mesh, connect, quant, time_step);

  /* The cell state is now updated at this stage. This will be useful for
     the monitoring */

  _update_binary_alloy_final_state(connect, quant, time_step);

  /* Update the forcing term in the momentum equation */

  alloy->update_velocity_forcing(mesh, connect, quant, time_step);

  if (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY)
    _compute_enthalpy(quant,
                      time_step->t_cur,        /* t_eval */
                      solid->temperature->val, /* temperature */
                      solid->g_l_field->val,   /* liquid fraction */
                      alloy->t_eut,            /* temp_ref */
                      solid->latent_heat,      /* latent heat coeff. */
                      solid->mass_density,     /* rho */
                      solid->cp,               /* cp */
                      solid->enthalpy->val);   /* computed enthalpy */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if solidification module is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_solidification_is_activated(void)
{
  if (cs_solidification_structure == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main structure to deal with solidification process
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_get_structure(void)
{
  return cs_solidification_structure;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the level of verbosity for the solidification module
 *
 * \param[in]   verbosity     level of verbosity to set
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_verbosity(int   verbosity)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  if (solid == NULL)
    return;

  solid->verbosity = verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value of the epsilon parameter used in the forcing term
 *         of the momentum equation
 *
 * \param[in]  forcing_eps    epsilon used in the penalization term to avoid a
 *                            division by zero
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_forcing_eps(cs_real_t    forcing_eps)
{
  assert(forcing_eps > 0);
  cs_solidification_forcing_eps = forcing_eps;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the solidification module
 *
 * \param[in]  model            type of modelling
 * \param[in]  options          flag to handle optional parameters
 * \param[in]  post_flag        predefined post-processings
 * \param[in]  boundaries       pointer to the domain boundaries
 * \param[in]  ns_model         model equations for the NavSto system
 * \param[in]  ns_model_flag    option flag for the Navier-Stokes system
 * \param[in]  algo_coupling    algorithm used for solving the NavSto system
 * \param[in]  ns_post_flag     predefined post-processings for Navier-Stokes
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_activate(cs_solidification_model_t       model,
                           cs_flag_t                       options,
                           cs_flag_t                       post_flag,
                           const cs_boundary_t            *boundaries,
                           cs_navsto_param_model_t         ns_model,
                           cs_navsto_param_model_flag_t    ns_model_flag,
                           cs_navsto_param_coupling_t      algo_coupling,
                           cs_navsto_param_post_flag_t     ns_post_flag)
{
  cs_flag_t  thm_num = 0, thm_post = 0, thm_model = 0;

  /* Allocate an empty structure */

  cs_solidification_t  *solid = _solidification_create();

  /* Set members of the structure according to the given settings */

  solid->model = model;

  /* By default, Stefan model is with a frozen field */

  if (model == CS_SOLIDIFICATION_MODEL_STEFAN)
    options |= CS_SOLIDIFICATION_NO_VELOCITY_FIELD;
  solid->options = options;

  if (post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS)
    post_flag |= CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE;
  solid->post_flag = post_flag;

  /* The Navier-Stokes is not solved when the frozen field is set */

  if (solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD) {

    solid->forcing_mom = NULL;
    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Incompatible set of options: "
                "no velocity and binary alloy model.\n"
                "Please check your settings.\n", __func__);

  }
  else {

    ns_model_flag |=
      CS_NAVSTO_MODEL_BOUSSINESQ | CS_NAVSTO_MODEL_WITH_SOLIDIFICATION;
    thm_model |= CS_THERMAL_MODEL_NAVSTO_ADVECTION;

    /* Add a property taking into account the head losses induced by the
       solidification process */

    solid->forcing_mom = cs_property_add("forcing_momentum_coef",
                                         CS_PROPERTY_ISO);

    /* If liquid, this coefficient is equal to zero */

    cs_property_set_reference_value(solid->forcing_mom, 0);

  }

  /* Activate the Navier-Stokes module */
  /* --------------------------------- */

  cs_navsto_system_t  *ns = cs_navsto_system_activate(boundaries,
                                                      ns_model,
                                                      ns_model_flag,
                                                      algo_coupling,
                                                      ns_post_flag);

  /* Activate and default settings for the thermal module */
  /* ---------------------------------------------------- */

  if (options & CS_SOLIDIFICATION_USE_ENTHALPY_VARIABLE)
    thm_model |= CS_THERMAL_MODEL_USE_ENTHALPY;
  else
    thm_model |= CS_THERMAL_MODEL_USE_TEMPERATURE;

  solid->thermal_sys = cs_thermal_system_activate(thm_model, thm_num, thm_post);

  cs_equation_param_t  *th_eqp =
    cs_equation_get_param(solid->thermal_sys->thermal_eq);

  /* Be sure that the space discretization is a CDO-Fb scheme */

  cs_equation_param_set(th_eqp, CS_EQKEY_SPACE_SCHEME, "cdofb");

  if (thm_model & CS_THERMAL_MODEL_USE_TEMPERATURE) {

    /* Add reaction property for the temperature equation */

    solid->thermal_reaction_coef = cs_property_add("thermal_reaction_coef",
                                                   CS_PROPERTY_ISO);

    cs_property_set_reference_value(solid->thermal_reaction_coef, 0);

    cs_equation_add_reaction(th_eqp, solid->thermal_reaction_coef);

  }

  /* Retrieve or add properties related to this module */

  solid->mass_density = ns->param->mass_density;
  assert(solid->mass_density != NULL);

  solid->viscosity = ns->param->tot_viscosity;
  assert(solid->viscosity != NULL);

  solid->cp = cs_property_by_name(CS_THERMAL_CP_NAME);
  assert(solid->cp != NULL);

  solid->g_l = cs_property_add("liquid_fraction", CS_PROPERTY_ISO);
  cs_property_set_reference_value(solid->g_l, 1.0);

  /* Initialize other members */

  solid->enthalpy = NULL;       /* Will be created later */
  solid->latent_heat = 1.0;     /* Default dummy value */

  /* Allocate the structure storing the modelling context/settings */

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_STEFAN:
    {
      cs_solidification_stefan_t  *s_model = NULL;
      BFT_MALLOC(s_model, 1, cs_solidification_stefan_t);

      /* Default initialization of this model */

      s_model->t_change = 0.;
      s_model->n_iter_max = 15;
      s_model->max_delta_h = 1e-2;

      /* Function pointers */

      solid->strategy = CS_SOLIDIFICATION_STRATEGY_PATH;
      s_model->update_gl = _update_gl_stefan;
      s_model->update_thm_st = _update_thm_stefan;

      /* Set the context */

      solid->model_context = (void *)s_model;

    }
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      cs_solidification_voller_t  *v_model = NULL;
      BFT_MALLOC(v_model, 1, cs_solidification_voller_t);

      /* Default initialization of this model */

      v_model->s_das = 0.33541;
      v_model->t_solidus = 0.;
      v_model->t_liquidus = 1.0;

      /* Non-linear algorithm */

      cs_iter_algo_param_t  nl_param = {
        .verbosity = 0,         /* level of display to output */
        .n_max_algo_iter = 15,  /* n_max iter. */
        .atol = 1e-6,           /* absolute tolerance */
        .rtol = 1e-2,           /* relative tolerance */
        .dtol = 1e3 };          /* divergence tolerance */

      if (solid->model == CS_SOLIDIFICATION_MODEL_VOLLER_NL) {

        v_model->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;
        v_model->nl_algo = cs_iter_algo_create(nl_param);

      }
      else {

        v_model->nl_algo_type = CS_PARAM_N_NL_ALGOS;
        v_model->nl_algo = NULL;

      }

      /* Function pointers */

      solid->strategy = CS_SOLIDIFICATION_STRATEGY_LEGACY;
      if (solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD)
        v_model->update_gl = _update_gl_voller_legacy_no_velocity;
      else
        v_model->update_gl = _update_gl_voller_legacy;
      v_model->update_thm_st = _update_thm_voller_legacy;

      /* If the non-linear model is used, then the default strategy is
         modified */

      if (solid->model == CS_SOLIDIFICATION_MODEL_VOLLER_NL) {
        solid->strategy = CS_SOLIDIFICATION_STRATEGY_PATH;
        v_model->update_thm_st = _update_thm_voller_path;
      }

      /* Set the context */

      solid->model_context = (void *)v_model;

    }
    break;

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      cs_solidification_binary_alloy_t *b_model = NULL;
      BFT_MALLOC(b_model, 1, cs_solidification_binary_alloy_t);

      /* Set a default value to the model parameters */

      b_model->ref_concentration = 1.;
      b_model->s_das = 0.33541;
      b_model->kp = 0.5;
      b_model->inv_kp = 2.;
      b_model->inv_kpm1 = -2;
      b_model->ml = -1;
      b_model->inv_ml = -1;
      b_model->t_melt = 0.;
      b_model->t_eut = b_model->t_eut_inf = b_model->t_eut_sup = 0.;
      b_model->c_eut = b_model->cs1 = b_model->dgldC_eut = 0.;

      b_model->diff_coef = 1e-6;

      /* Monitoring and criteria to drive the convergence of the coupled system
         (solute transport and thermal equation) */

      b_model->iter = 0;
      b_model->n_iter_max = 10;
      b_model->delta_tolerance = 1e-3;
      b_model->eta_relax = 0.;
      b_model->gliq_relax = 0.;

      /* Strategy to perform the main steps of the simulation of a binary alloy
       * Default strategy: Taylor which corresponds to the Legacy one with
       * improvements thanks to some Taylor expansions */

      solid->strategy = CS_SOLIDIFICATION_STRATEGY_TAYLOR;

      /* Functions which are specific to a strategy */

      b_model->update_gl = _update_gl_taylor;
      b_model->update_thm_st = _update_thm_taylor;

      /* Functions which are common to all strategies */

      b_model->update_velocity_forcing = _update_velocity_forcing;
      b_model->update_clc = _update_clc;

      /* Initialize pointers */

      b_model->tk_bulk = NULL;
      b_model->ck_bulk = NULL;
      b_model->tx_bulk = NULL;
      b_model->cx_bulk = NULL;

      /* alloy->temp_faces is shared and set when defined */

      b_model->c_l_cells = NULL;
      b_model->c_l_faces = NULL;

      b_model->solute_equation = NULL;
      b_model->c_bulk = NULL;

      b_model->diff_pty = NULL;
      b_model->diff_pty_array = NULL;

      /* eta_coef is activated with advanced options */

      b_model->eta_coef_pty = NULL;
      b_model->eta_coef_array = NULL;

      /* Optional post-processing arrays */

      b_model->t_liquidus = NULL;
      b_model->tbulk_minus_tliq = NULL;
      b_model->cliq_minus_cbulk = NULL;

      /* Set the context */

      solid->model_context = (void *)b_model;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model for the solidification module.\n"
              " Please check your setup.", __func__);

  } /* Switch on the solidification model */

  /* Set the global pointer */

  cs_solidification_structure = solid;

  return solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the Stefan model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_stefan_t *
cs_solidification_get_stefan_struct(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  if (solid->model != CS_SOLIDIFICATION_MODEL_STEFAN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Stefan model not declared during the"
              " activation of the solidification module.\n"
              " Please check your settings.", __func__);

  cs_solidification_stefan_t  *s_model
    = (cs_solidification_stefan_t *)solid->model_context;
  assert(s_model != NULL);

  return s_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the Stefan's model settings
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_stefan_t *
cs_solidification_check_stefan_model(void)
{
  cs_solidification_stefan_t  *s_model = cs_solidification_get_stefan_struct();

  if (s_model->n_iter_max < 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n_iter_max (= %d)\n",
              __func__, s_model->n_iter_max);
  if (s_model->max_delta_h < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for max_delta_h (= %5.3e)\n",
              __func__, s_model->max_delta_h);

  return s_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which describe the Stefan model
 *
 * \param[in] t_change     liquidus/solidus temperature (in K)
 * \param[in] latent_heat  latent heat
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_stefan_model(cs_real_t    t_change,
                                   cs_real_t    latent_heat)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_stefan_t  *s_model = cs_solidification_get_stefan_struct();

  /* Model parameters */

  s_model->t_change = t_change;
  solid->latent_heat = latent_heat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the Voller model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_voller_t *
cs_solidification_get_voller_struct(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  if (solid->model != CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 &&
      solid->model != CS_SOLIDIFICATION_MODEL_VOLLER_NL )
    bft_error(__FILE__, __LINE__, 0,
              " %s: Voller model not declared during the"
              " activation of the solidification module.\n"
              " Please check your settings.", __func__);

  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;
  assert(v_model != NULL);

  return v_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the Voller's model settings
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_voller_t *
cs_solidification_check_voller_model(void)
{
  cs_solidification_voller_t  *v_model = cs_solidification_get_voller_struct();

  if (v_model->t_liquidus - v_model->t_solidus < 0.)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The liquidus and solidus temperatures are not"
              " consistent.\n"
              " Please check your settings.", __func__);
  if (v_model->s_das < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for the secondary dendrite arms spacing",
              __func__, v_model->s_das);

  return v_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which describe the Voller and
 *         Prakash modelling
 *
 * \param[in]  beta           thermal dilatation coefficient
 * \param[in]  t_ref          reference temperature (for the Boussinesq approx)
 * \param[in]  t_solidus      solidus temperature (in K)
 * \param[in]  t_liquidus     liquidus temperature (in K)
 * \param[in]  latent_heat    latent heat
 * \param[in]  s_das          secondary dendrite space arms
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_voller_model(cs_real_t    beta,
                                   cs_real_t    t_ref,
                                   cs_real_t    t_solidus,
                                   cs_real_t    t_liquidus,
                                   cs_real_t    latent_heat,
                                   cs_real_t    s_das)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model = cs_solidification_get_voller_struct();

  /* Model parameters */

  v_model->t_solidus = t_solidus;
  v_model->t_liquidus = t_liquidus;
  v_model->s_das = s_das;

  solid->latent_heat = latent_heat;

  /* Add the Boussinesq term */

  cs_navsto_param_t *nsp = cs_navsto_system_get_param();

  cs_navsto_param_add_boussinesq_term(nsp, beta, t_ref);

  /* Set the reference temperature in the thermal module */

  cs_thermal_system_set_reference_temperature(t_ref);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which describe the Voller and
 *         Prakash modelling
 *
 * \param[in]  beta           thermal dilatation coefficient
 * \param[in]  t_ref          reference temperature (for the Boussinesq approx)
 * \param[in]  t_solidus      solidus temperature (in K)
 * \param[in]  t_liquidus     liquidus temperature (in K)
 * \param[in]  latent_heat    latent heat
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_voller_model_no_velocity(cs_real_t    t_solidus,
                                               cs_real_t    t_liquidus,
                                               cs_real_t    latent_heat)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model = cs_solidification_get_voller_struct();

  if ((solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD) == 0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: CS_SOLIDIFICATION_NO_VELOCITY_FIELD has not been set with"
              " the Voller model.\n"
              "Please check your settings.\n", __func__);

  /* Model parameters (those which are useful in this case) */

  v_model->t_solidus = t_solidus;
  v_model->t_liquidus = t_liquidus;
  solid->latent_heat = latent_heat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the binary alloy model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_binary_alloy_t *
cs_solidification_get_binary_alloy_struct(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  if (solid->model != CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The binary alloy model was not declared during the"
              " activation of the solidification module.\n"
              " Please check your settings.", __func__);

  cs_solidification_binary_alloy_t  *b_model
    = (cs_solidification_binary_alloy_t *)solid->model_context;
  assert(b_model != NULL);

  return b_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the settings of the binary alloy
 *         model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_binary_alloy_t *
cs_solidification_check_binary_alloy_model(void)
{
  cs_solidification_binary_alloy_t
    *b_model = cs_solidification_get_binary_alloy_struct();

  if (b_model->s_das < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for the secondary dendrite arms spacing",
              __func__, b_model->s_das);
  if (b_model->kp < FLT_MIN || b_model->kp > 1 - FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for partition coefficient",
              __func__, b_model->kp);
  if (fabs(b_model->ml) < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for the liquidus slope",
              __func__, b_model->ml);
  if (b_model->n_iter_max < 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for n_iter_max (current: %d).\n"
              " Should be strictly greater than 0.\n",
              __func__, b_model->n_iter_max);
  if (b_model->delta_tolerance < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid value for \"tolerance\" (current: %5.3e).\n",
              __func__, b_model->delta_tolerance);

  return b_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which describe a solidification
 *         process with a binary alloy (with components A and B)
 *         Add a transport equation for the solute concentration to simulate
 *         the conv/diffusion of the alloy ratio between the two components of
 *         the alloy
 *
 * \param[in]  name          name of the binary alloy
 * \param[in]  varname       name of the unknown related to the tracer eq.
 * \param[in]  beta_t        thermal dilatation coefficient
 * \param[in]  temp0         reference temperature (Boussinesq term)
 * \param[in]  beta_c        solutal dilatation coefficient
 * \param[in]  conc0         reference mixture concentration (Boussinesq term)
 * \param[in]  kp            value of the distribution coefficient
 * \param[in]  mliq          liquidus slope for the solute concentration
 * \param[in]  t_eutec       temperature at the eutectic point
 * \param[in]  t_melt        phase-change temperature for the pure material (A)
 * \param[in]  solute_diff   solutal diffusion coefficient in the liquid
 * \param[in]  latent_heat   latent heat
 * \param[in]  s_das         secondary dendrite arm spacing
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_binary_alloy_model(const char     *name,
                                         const char     *varname,
                                         cs_real_t       beta_t,
                                         cs_real_t       temp0,
                                         cs_real_t       beta_c,
                                         cs_real_t       conc0,
                                         cs_real_t       kp,
                                         cs_real_t       mliq,
                                         cs_real_t       t_eutec,
                                         cs_real_t       t_melt,
                                         cs_real_t       solute_diff,
                                         cs_real_t       latent_heat,
                                         cs_real_t       s_das)
{
  /* Check the validity of some parameters */

  if (kp < FLT_MIN || kp > 1 - FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for partition coefficient", __func__, kp);
  if (fabs(mliq) < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for the liquidus slope", __func__, mliq);
  if (s_das < FLT_MIN)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid value %g for the secondary dendrite arms spacing",
              __func__, s_das);

  /* Retrieve and set the binary alloy structures */

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_binary_alloy_t
    *alloy = cs_solidification_get_binary_alloy_struct();

  solid->latent_heat = latent_heat;

  /* Set the main physical parameters/constants */

  alloy->ref_concentration = conc0;
  alloy->s_das = s_das;

  /* Phase diagram parameters and related quantities */

  alloy->kp = kp;               /* partition coeff. */
  alloy->ml = mliq;             /* liquidus slope */
  alloy->t_melt = t_melt;       /* melting temperature */
  alloy->t_eut = t_eutec;       /* eutectic temperature */

  /* Derived quantities from the phase diagram */

  alloy->inv_kp = 1./kp;
  alloy->inv_kpm1 = 1./(alloy->kp - 1.);
  alloy->inv_ml = 1./mliq;
  alloy->c_eut = (t_eutec - t_melt)*alloy->inv_ml;
  alloy->cs1 = alloy->c_eut * kp; /* Apply the lever rule */

  assert(fabs(alloy->c_eut - alloy->cs1) > FLT_MIN);
  alloy->dgldC_eut = 1./(alloy->c_eut - alloy->cs1);

  /* Define a small range of temperature around the eutectic temperature in
   * which one assumes an eutectic transformation */

  alloy->t_eut_inf =
    alloy->t_eut - cs_solidification_eutectic_threshold;
  alloy->t_eut_sup =
    alloy->t_eut + cs_solidification_eutectic_threshold;

  /* Alloy equation and variable field */

  assert(name != NULL && varname != NULL);
  cs_equation_t  *eq = cs_equation_add(name, varname,
                                       CS_EQUATION_TYPE_SOLIDIFICATION,
                                       1,
                                       CS_PARAM_BC_HMG_NEUMANN);

  /* Set an upwind scheme by default since it could be a pure advection eq. */

  cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  /* Set the default numerical options that should be used */

  cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
  cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
  cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "sushi");
  cs_equation_param_set(eqp, CS_EQKEY_ADV_SCHEME, "upwind");
  cs_equation_param_set(eqp, CS_EQKEY_ADV_FORMULATION, "conservative");

  alloy->solute_equation = eq;
  alloy->c_bulk = NULL;  /* Variable field related to this equation. This will
                            be set later (after the equation initialization) */

  /* Always add a diffusion term (to avoid a zero block face-face when there
     is no more convection */

  if (solute_diff > 0)
    alloy->diff_coef = solute_diff;
  else
    alloy->diff_coef = cs_solidification_diffusion_eps;

  char  *pty_name = NULL;
  size_t  len = strlen(varname) + strlen("_diff_pty");
  BFT_MALLOC(pty_name, len + 1, char);
  sprintf(pty_name, "%s_diff_pty", varname);
  pty_name[len] = '\0';
  alloy->diff_pty = cs_property_add(pty_name, CS_PROPERTY_ISO);
  BFT_FREE(pty_name);

  cs_equation_add_diffusion(eqp, alloy->diff_pty);

  /* Add Boussinesq terms (two parts) */

  cs_navsto_param_t *nsp = cs_navsto_system_get_param();

  /* Thermal effect */

  cs_navsto_param_add_boussinesq_term(nsp, beta_t, temp0);

  /* Solutal effect */

  cs_navsto_param_add_boussinesq_term(nsp, beta_c, conc0);

  /* Set the reference temperature in the thermal module */

  cs_thermal_system_set_reference_temperature(temp0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the strategy to update quantitiess (liquid fraction and
 *         the thermal source term for the two main quantities)
 *
 * \param[in]  strategy     strategy to perform the update of quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_strategy(cs_solidification_strategy_t  strategy)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_STEFAN:
    cs_base_warn(__FILE__, __LINE__);
    bft_printf("%s:  Only one strategy is available with the Stefan model.\n",
               __func__);
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      cs_solidification_voller_t  *v_model =
        cs_solidification_get_voller_struct();

      switch (strategy) {

      case CS_SOLIDIFICATION_STRATEGY_LEGACY:
        v_model->update_thm_st = _update_thm_voller_legacy;
        break;

      case CS_SOLIDIFICATION_STRATEGY_PATH:
        v_model->update_thm_st = _update_thm_voller_path;
        break;

      default:
        if (solid->model == CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87)
          v_model->update_thm_st = _update_thm_voller_legacy;
        else
          v_model->update_thm_st = _update_thm_voller_path;
        break;

      } /* Switch on the strategy */

    }
    break; /* Voller-like models */

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      cs_solidification_binary_alloy_t *alloy =
        cs_solidification_get_binary_alloy_struct();

      switch (strategy) {

      case CS_SOLIDIFICATION_STRATEGY_LEGACY:
        if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM)
          alloy->update_gl = _update_gl_legacy_ast;
        else
          alloy->update_gl = _update_gl_legacy;
        alloy->update_thm_st = _update_thm_legacy;
        break;

      case CS_SOLIDIFICATION_STRATEGY_TAYLOR:
        if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Adding an advective source term is incompatible with"
                    " the Taylor strategy.\n", __func__);
        else
          alloy->update_gl = _update_gl_taylor;
        alloy->update_thm_st = _update_thm_taylor;
        break;

      case CS_SOLIDIFICATION_STRATEGY_PATH:
        if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Adding an advective source term is incompatible with"
                    " the Path strategy.\n", __func__);
        else
          alloy->update_gl = _update_gl_binary_path;
        alloy->update_thm_st = _update_thm_binary_path;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, "%s: Invalid strategy.\n", __func__);
        break;

      } /* Switch on strategies */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid solidification model.\n", __func__);

  } /* Switch on the solidification model */

  solid->strategy = strategy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the functions to perform the update of physical properties
 *         and/or the computation of the thermal source term or quantities
 *         and/or the way to perform the coupling between the thermal equation
 *         and the bulk concentration computation. All this setting defines
 *         the way to compute the solidification process of a binary alloy.
 *         If a function is set to NULL then the automatic settings are kept.
 *
 *         --Advanced usage-- This enables to finely control the numerical or
 *         physical modelling aspects.
 *
 * \param[in] vel_forcing        pointer to update the velocity forcing
 * \param[in] cliq_update        pointer to update the liquid concentration
 * \param[in] gliq_update        pointer to update the liquid fraction
 * \param[in] thm_st_update      pointer to update thermal source terms
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_segr_functions(cs_solidification_func_t  *vel_forcing,
                                     cs_solidification_func_t  *cliq_update,
                                     cs_solidification_func_t  *gliq_update,
                                     cs_solidification_func_t  *thm_st_update)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  cs_solidification_binary_alloy_t  *alloy
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  assert(solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY);
  assert(alloy != NULL);

  if (vel_forcing != NULL) {
    alloy->update_velocity_forcing = vel_forcing;
    solid->options |= CS_SOLIDIFICATION_BINARY_ALLOY_M_FUNC;
  }

  if (cliq_update != NULL) {
    alloy->update_clc = cliq_update;
    solid->options |= CS_SOLIDIFICATION_BINARY_ALLOY_C_FUNC;
  }

  if (gliq_update != NULL) {
    alloy->update_gl = gliq_update;
    solid->options |= CS_SOLIDIFICATION_BINARY_ALLOY_G_FUNC;
  }

  if (thm_st_update != NULL) {
    alloy->update_thm_st = thm_st_update;
    solid->options |= CS_SOLIDIFICATION_BINARY_ALLOY_T_FUNC;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the solidification module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_destroy_all(void)
{
  if (cs_solidification_structure == NULL)
    return NULL;

  cs_solidification_t  *solid = cs_solidification_structure;

  /* The lifecycle of properties, equations and fields is not managed by
   * the current structure and sub-structures.
   * Free only what is owned by this structure */

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_STEFAN:
    {
      cs_solidification_stefan_t  *s_model
        = (cs_solidification_stefan_t *)solid->model_context;

      BFT_FREE(s_model);

    } /* Stefan modelling */
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      cs_solidification_voller_t  *v_model
        = (cs_solidification_voller_t *)solid->model_context;

      if (CS_SOLIDIFICATION_MODEL_VOLLER_NL)
        BFT_FREE(v_model->nl_algo);

      BFT_FREE(v_model);

    } /* Voller and Prakash modelling */
    break;

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      cs_solidification_binary_alloy_t  *alloy
        = (cs_solidification_binary_alloy_t *)solid->model_context;

      BFT_FREE(alloy->diff_pty_array);
      BFT_FREE(alloy->c_l_cells);
      BFT_FREE(alloy->eta_coef_array);
      BFT_FREE(alloy->tk_bulk);
      BFT_FREE(alloy->ck_bulk);

      if (solid->options & CS_SOLIDIFICATION_USE_EXTRAPOLATION) {
        BFT_FREE(alloy->tx_bulk);
        BFT_FREE(alloy->cx_bulk);
      }

      if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM)
        BFT_FREE(alloy->c_l_faces);

      if (solid->post_flag & CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE)
        BFT_FREE(alloy->t_liquidus);

      if (solid->post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS) {
        BFT_FREE(alloy->tbulk_minus_tliq);
        BFT_FREE(alloy->cliq_minus_cbulk);
      }

      BFT_FREE(alloy);

    } /* Binary alloy modelling */
    break;

  default:
    break; /* Nothing to do */

  } /* Switch on solidification model */

  BFT_FREE(solid->thermal_reaction_coef_array);
  BFT_FREE(solid->thermal_source_term_array);
  BFT_FREE(solid->forcing_mom_array);

  BFT_FREE(solid->cell_state);

  if (solid->plot_state != NULL)
    cs_time_plot_finalize(&solid->plot_state);

  BFT_FREE(solid);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the solidification module
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_setup(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");

  /* Add a field for the liquid fraction */

  solid->g_l_field = cs_field_create("liquid_fraction",
                                     field_mask,
                                     c_loc_id,
                                     1,
                                     true); /* has_previous */

  cs_field_set_key_int(solid->g_l_field, log_key, 1);
  cs_field_set_key_int(solid->g_l_field, post_key, 1);

  /* Add the enthalpy field if not already created */

  solid->enthalpy = cs_field_by_name_try("enthalpy");
  if (solid->enthalpy == NULL) {

    bool add_enthalpy = false;

    if (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY)
      add_enthalpy = true;
    if (solid->model == CS_SOLIDIFICATION_MODEL_VOLLER_NL ||
        solid->model == CS_SOLIDIFICATION_MODEL_STEFAN)
      add_enthalpy = true;

    if (add_enthalpy) {
      solid->enthalpy = cs_field_create("enthalpy",
                                        field_mask,
                                        c_loc_id,
                                        1,
                                        true); /* has_previous */

      cs_field_set_key_int(solid->enthalpy, log_key, 1);
    }

  }

  if (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY)
    cs_field_set_key_int(solid->enthalpy, post_key, 1);

  /* Add a reaction term to the momentum equation */

  if (solid->forcing_mom != NULL) {

    cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
    cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
    assert(mom_eqp != NULL);

    cs_equation_add_reaction(mom_eqp, solid->forcing_mom);

  }

  /* Add default post-processing related to the solidification module */

  cs_post_add_time_mesh_dep_output(cs_solidification_extra_post, solid);

  /* Model-specific part */

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      /* Check the sanity of the model parameters and retrieve the structure */

      cs_solidification_voller_t
        *v_model = cs_solidification_check_voller_model();

      solid->forcing_coef = 180./(v_model->s_das*v_model->s_das);
    }
    break;

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      /* Check the sanity of the model parameters and retrieve the structure */

      cs_solidification_binary_alloy_t
        *alloy = cs_solidification_check_binary_alloy_model();

      cs_equation_param_t  *eqp = cs_equation_get_param(alloy->solute_equation);

      /* Add the unsteady term */

      cs_equation_add_time(eqp, solid->mass_density);

      /* Add an advection term to the solute concentration equation */

      cs_equation_add_advection(eqp, cs_navsto_get_adv_field());

      if ((solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM) == 0) {

        alloy->eta_coef_pty = cs_property_add("alloy_adv_coef",
                                              CS_PROPERTY_ISO);

        cs_equation_add_advection_scaling_property(eqp, alloy->eta_coef_pty);

      }

      solid->forcing_coef = 180./(alloy->s_das*alloy->s_das);

      /* Add the variable field (automatic) */

      cs_equation_predefined_create_field(-1, alloy->solute_equation);
    }
    break; /* Binary alloy model */

  default: /* Stefan: There is nothing else to do */
    break;

  } /* Switch on model */

  if (cs_glob_rank_id < 1) {

    int  n_output_states = CS_SOLIDIFICATION_N_STATES - 1;
    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
      n_output_states += 1;

    int  n_output_values = n_output_states;
    if (solid->post_flag & CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE)
      n_output_values += 1;

    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {
      if (solid->post_flag & CS_SOLIDIFICATION_POST_SEGREGATION_INDEX)
        n_output_values += 1;
    }

    const char  **labels;
    BFT_MALLOC(labels, n_output_values, const char *);
    for (int i = 0; i < n_output_states; i++)
      labels[i] = _state_names[i];

    n_output_values = n_output_states;
    if (solid->post_flag & CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE)
      labels[n_output_values++] = "SolidRate";

    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY)
      if (solid->post_flag & CS_SOLIDIFICATION_POST_SEGREGATION_INDEX)
        labels[n_output_values++] = "SegrIndex";

    /* Use the physical time rather than the number of iterations */

    if (n_output_values > 0)
      solid->plot_state = cs_time_plot_init_probe("solidification",
                                                  "",
                                                  CS_TIME_PLOT_DAT,
                                                  false,
                                                  180,   /* flush time */
                                                  -1,
                                                  n_output_values,
                                                  NULL,
                                                  NULL,
                                                  labels);

    BFT_FREE(labels);

  } /* rank 0 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the solidification
 *         module
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_finalize_setup(const cs_cdo_connect_t       *connect,
                                 const cs_cdo_quantities_t    *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  const cs_lnum_t  n_cells = quant->n_cells;
  const size_t  size_c = n_cells*sizeof(cs_real_t);

  /* Retrieve the field associated to the temperature */

  solid->temperature = cs_field_by_name("temperature");

  /* Define the liquid fraction */

  cs_property_def_by_field(solid->g_l, solid->g_l_field);

  /* Initially one assumes that all is liquid except for cells in a
   * predefined solid zone for all the computation */

  BFT_MALLOC(solid->cell_state, n_cells, cs_solidification_state_t);

  cs_field_set_values(solid->g_l_field, 1.);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++) {

    if (connect->cell_flag[i] & CS_FLAG_SOLID_CELL) {
      solid->g_l_field->val[i] = 0;
      solid->g_l_field->val_pre[i] = 0;
      solid->cell_state[i] = CS_SOLIDIFICATION_STATE_SOLID;
    }
    else {
      solid->g_l_field->val_pre[i] = 1.;
      solid->cell_state[i] = CS_SOLIDIFICATION_STATE_LIQUID;
    }

  } /* Loop on cells */

  if (cs_flag_test(solid->options,
                   CS_SOLIDIFICATION_NO_VELOCITY_FIELD) == false) {

    /* Define the forcing term acting as a reaction term in the momentum
       equation. This term is related to the liquid fraction */

    BFT_MALLOC(solid->forcing_mom_array, n_cells, cs_real_t);
    memset(solid->forcing_mom_array, 0, size_c);

    cs_property_def_by_array(solid->forcing_mom,
                             cs_flag_primal_cell,
                             solid->forcing_mom_array,
                             false, /* definition is owner ? */
                             NULL, NULL); /* no index, no ids */

    /* Add the temperature array for the Boussinesq term (thermal effect) */

    cs_navsto_param_t *nsp = cs_navsto_system_get_param();

    assert(nsp->n_boussinesq_terms > 0);
    cs_navsto_param_boussinesq_t  *bp = nsp->boussinesq_param;

    cs_navsto_param_set_boussinesq_array(bp, solid->temperature->val);

  }

  /* Define the reaction coefficient and the source term for the temperature
     equation */

  if (solid->thermal_reaction_coef != NULL) {

    BFT_MALLOC(solid->thermal_reaction_coef_array, n_cells, cs_real_t);
    memset(solid->thermal_reaction_coef_array, 0, size_c);

    cs_property_def_by_array(solid->thermal_reaction_coef,
                             cs_flag_primal_cell,
                             solid->thermal_reaction_coef_array,
                             false, /* definition is owner ? */
                             NULL, NULL); /* no index, no ids */

    BFT_MALLOC(solid->thermal_source_term_array, n_cells, cs_real_t);
    memset(solid->thermal_source_term_array, 0, size_c);

    cs_equation_param_t  *thm_eqp =
      cs_equation_param_by_name(CS_THERMAL_EQNAME);

    cs_equation_add_source_term_by_array(thm_eqp,
                                         NULL,   /* all cells selected */
                                         cs_flag_primal_cell,
                                         solid->thermal_source_term_array,
                                         false,  /* definition is owner ? */
                                         NULL, NULL); /* no index, no ids */

  }

  if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {
    /*                ==================================== */

    cs_solidification_binary_alloy_t  *alloy
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    /* Get a shortcut to the c_bulk field */

    alloy->c_bulk = cs_equation_get_field(alloy->solute_equation);

    /* Allocate an array to store the liquid concentration */

    BFT_MALLOC(alloy->c_l_cells, n_cells, cs_real_t);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      alloy->c_l_cells[i] = alloy->ref_concentration;

    /* Add the c_l_cells array for the Boussinesq term (solutal effect) */

    cs_navsto_param_t *nsp = cs_navsto_system_get_param();

    assert(nsp->n_boussinesq_terms == 2);
    cs_navsto_param_boussinesq_t  *bp = nsp->boussinesq_param + 1;

    cs_navsto_param_set_boussinesq_array(bp, alloy->c_l_cells);

    /* Allocate arrays storing the intermediate states during the
       sub-iterations to solve the non-linearity */

    BFT_MALLOC(alloy->tk_bulk, n_cells, cs_real_t);
    BFT_MALLOC(alloy->ck_bulk, n_cells, cs_real_t);

    if (solid->options & CS_SOLIDIFICATION_USE_EXTRAPOLATION) {
      BFT_MALLOC(alloy->tx_bulk, n_cells, cs_real_t);
      BFT_MALLOC(alloy->cx_bulk, n_cells, cs_real_t);
    }

    /* Allocate eta even if SOLUTE_WITH_SOURCE_TERM is activated */

    const cs_real_t  eta_ref_value = 1.;
    BFT_MALLOC(alloy->eta_coef_array, n_cells, cs_real_t);
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      alloy->eta_coef_array[i] = eta_ref_value;

    if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM) {

      BFT_MALLOC(alloy->c_l_faces, quant->n_faces, cs_real_t);
      memset(alloy->c_l_faces, 0, sizeof(cs_real_t)*quant->n_faces);

    }
    else { /* Estimate the reference value for the solutal diffusion property
            * One assumes that g_l (the liquid fraction is equal to 1) */

      cs_property_set_reference_value(alloy->eta_coef_pty, eta_ref_value);

      cs_property_def_by_array(alloy->eta_coef_pty,
                               cs_flag_primal_cell,
                               alloy->eta_coef_array,
                               false,
                               NULL, NULL); /* no index, no ids */

    }

    /* Estimate the reference value for the solutal diffusion property
     * One assumes that g_l (the liquid fraction) is equal to 1. */

    const cs_real_t  pty_ref_value =
      solid->mass_density->ref_value*alloy->diff_coef;

    cs_property_set_reference_value(alloy->diff_pty, pty_ref_value);

    BFT_MALLOC(alloy->diff_pty_array, n_cells, cs_real_t);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      alloy->diff_pty_array[i] = pty_ref_value;

    cs_property_def_by_array(alloy->diff_pty,
                             cs_flag_primal_cell,
                             alloy->diff_pty_array,
                             false,
                             NULL, NULL); /* no index/ids */

    if (solid->post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS) {

      BFT_MALLOC(alloy->tbulk_minus_tliq, n_cells, cs_real_t);
      memset(alloy->tbulk_minus_tliq, 0, size_c);
      BFT_MALLOC(alloy->cliq_minus_cbulk, n_cells, cs_real_t);
      memset(alloy->cliq_minus_cbulk, 0, size_c);

    }

    if (solid->post_flag & CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE)
      BFT_MALLOC(alloy->t_liquidus, n_cells, cs_real_t);

  } /* Binary alloy model */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the solidification module in the log file dedicated to
 *         the setup
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_log_setup(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL)
    return;

  const char  *module = "Solidification";

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the solidification module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  cs_log_printf(CS_LOG_SETUP, "  * %s | Verbosity: %d\n",
                module, solid->verbosity);

  /* Display options */

  cs_log_printf(CS_LOG_SETUP, "  * %s | Strategy:", module);
  switch (solid->strategy) {

  case CS_SOLIDIFICATION_STRATEGY_LEGACY:
    cs_log_printf(CS_LOG_SETUP, " **Legacy**\n");
    break;
  case CS_SOLIDIFICATION_STRATEGY_TAYLOR:
    cs_log_printf(CS_LOG_SETUP, " **Legacy + Taylor-based updates**\n");
    break;
  case CS_SOLIDIFICATION_STRATEGY_PATH:
    cs_log_printf(CS_LOG_SETUP, " **Rely on the solidification path**\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid strategy\n", __func__);
  }

  switch (solid->model) {
  case CS_SOLIDIFICATION_MODEL_STEFAN:
    {
      cs_solidification_stefan_t  *s_model =
        (cs_solidification_stefan_t *)solid->model_context;

      cs_log_printf(CS_LOG_SETUP, "  * %s | Model: **Stefan**\n", module);
      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Tliq/sol: %5.3e\n"
                    "  * %s | Latent heat: %5.3e\n"
                    "  * %s | Max. iter: %d; Max. delta enthalpy: %5.3e\n",
                    module, s_model->t_change, module, solid->latent_heat,
                    module, s_model->n_iter_max, s_model->max_delta_h);
    }
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
    {
      cs_solidification_voller_t  *v_model
        = (cs_solidification_voller_t *)solid->model_context;

      cs_log_printf(CS_LOG_SETUP, "  * %s | Model: **Voller-Prakash (1987)**\n",
                    module);
      if (solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD)
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Tliq: %5.3e; Tsol: %5.3e\n"
                      "  * %s | Latent heat: %5.3e\n",
                      module, v_model->t_liquidus, v_model->t_solidus,
                      module, solid->latent_heat);
      else
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Tliq: %5.3e; Tsol: %5.3e\n"
                      "  * %s | Latent heat: %5.3e\n"
                      "  * %s | Forcing coef: %5.3e s_das: %5.3e\n",
                      module, v_model->t_liquidus, v_model->t_solidus,
                      module, solid->latent_heat,
                      module, solid->forcing_coef, v_model->s_das);

    }
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      cs_solidification_voller_t  *v_model
        = (cs_solidification_voller_t *)solid->model_context;

      cs_iter_algo_t  *ia = v_model->nl_algo;
      assert(ia != NULL);

      cs_log_printf(CS_LOG_SETUP, "  * %s |"
                    " **Model: Voller-Prakash (1987) with non-linearities**\n",
                    module);

      if (solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD)
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Tliq: %5.3e; Tsol: %5.3e\n"
                      "  * %s | Latent heat: %5.3e\n"
                      "  * %s | NL Algo: max. iter: %d; rtol: %5.3e,"
                      " atol: %5.3e, dtol: %5.3e\n",
                      module, v_model->t_liquidus, v_model->t_solidus,
                      module, solid->latent_heat,
                      module, ia->param.n_max_algo_iter, ia->param.rtol,
                      ia->param.atol, ia->param.dtol);
      else
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Tliq: %5.3e; Tsol: %5.3e\n"
                      "  * %s | Latent heat: %5.3e\n"
                      "  * %s | Forcing coef: %5.3e s_das: %5.3e\n"
                      "  * %s | NL Algo: max. iter: %d; rtol: %5.3e,"
                      " atol: %5.3e, dtol: %5.3e\n",
                      module, v_model->t_liquidus, v_model->t_solidus,
                      module, solid->latent_heat,
                      module, solid->forcing_coef, v_model->s_das,
                      module, ia->param.n_max_algo_iter, ia->param.rtol,
                      ia->param.atol, ia->param.dtol);
    }
    break;

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      cs_solidification_binary_alloy_t  *alloy
        = (cs_solidification_binary_alloy_t *)solid->model_context;

      cs_log_printf(CS_LOG_SETUP, "  * %s | Model: **Binary alloy**\n",
                    module);
      cs_log_printf(CS_LOG_SETUP, "  * %s | Alloy: %s\n",
                    module, cs_equation_get_name(alloy->solute_equation));
      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Distribution coef.: %5.3e\n"
                    "  * %s | Liquidus slope: %5.3e\n"
                    "  * %s | Phase change temp.: %5.3e\n"
                    "  * %s | Eutectic conc.: %5.3e\n"
                    "  * %s | Latent heat: %5.3e\n",
                    module, alloy->kp,
                    module, alloy->ml, module, alloy->t_melt,
                    module, alloy->c_eut,
                    module, solid->latent_heat);
      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Forcing coef: %5.3e; s_das: %5.3e\n",
                    module, solid->forcing_coef, alloy->s_das);

      cs_log_printf(CS_LOG_SETUP, "  * %s | Options:", module);
      if (solid->options & CS_SOLIDIFICATION_BINARY_ALLOY_C_FUNC)
        cs_log_printf(CS_LOG_SETUP,
                      " User-defined function for the concentration eq.");
      else {

        if (cs_flag_test(solid->options,
                         CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM))
          cs_log_printf(CS_LOG_SETUP,
                        " Solute concentration with an advective source term");
        else
          cs_log_printf(CS_LOG_SETUP,
                        " Solute concentration with an advective coefficient");

      } /* Not user-defined */
      cs_log_printf(CS_LOG_SETUP, "\n");

      if (solid->options & CS_SOLIDIFICATION_BINARY_ALLOY_T_FUNC)
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Options: %s\n", module,
                      " User-defined function for the thermal equation");

      if (solid->options & CS_SOLIDIFICATION_BINARY_ALLOY_G_FUNC)
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Options: %s\n", module,
                      " User-defined function for the liquid fraction/state");

      if (cs_flag_test(solid->options, CS_SOLIDIFICATION_USE_EXTRAPOLATION))
        cs_log_printf(CS_LOG_SETUP,
                      "  * %s | Options: %s\n", module,
                      " Update using a second-order in time extrapolation");

      if (solid->options & CS_SOLIDIFICATION_WITH_PENALIZED_EUTECTIC) {
        if (solid->strategy == CS_SOLIDIFICATION_STRATEGY_PATH)
          cs_log_printf(CS_LOG_SETUP, "  * %s | Options: %s\n", module,
                      " Penalized eutectic temperature");
      else
        cs_log_printf(CS_LOG_SETUP, "  * %s | Options: %s\n", module,
                      " Penalized eutectic temperature (unused)");
      }

      if (alloy->n_iter_max > 1)
        cs_log_printf(CS_LOG_SETUP, "  * %s | Options: Use sub-iterations"
                      " n_iter_max %d; tolerance: %.3e\n",
                      module, alloy->n_iter_max, alloy->delta_tolerance);

    } /* Binary alloy */

  default:
    break;

  } /* Switch on model type */

  cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set an initial values for all quantities related to this module
 *         This is done after the setup step.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_values(const cs_mesh_t              *mesh,
                              const cs_cdo_connect_t       *connect,
                              const cs_cdo_quantities_t    *quant,
                              const cs_time_step_t         *time_step)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);

  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  /* Set the first fluid/solid cell and sanity check for the mass density in the
     fluid/solid zone */

  const cs_real_t  cp0 = solid->cp->ref_value;
  const cs_real_t  rho0 = solid->mass_density->ref_value;

  for (int i = 0; i < cs_volume_zone_n_zones(); i++) {

    const cs_zone_t  *z = cs_volume_zone_by_id(i);

    if (z->type & CS_VOLUME_ZONE_SOLID) /* permanent solid zone */
      continue;

    else { /* fluid/solid zone according to thermodynamics conditions */

      if (z->n_elts == 0)
        continue;

      if (solid->first_cell < 0)
        solid->first_cell = z->elt_ids[0];

      else {

        cs_real_t  rho = cs_property_get_cell_value(z->elt_ids[0],
                                                    time_step->t_cur,
                                                    solid->mass_density);

        if (fabs(rho - rho0) > FLT_MIN)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: A uniform value of the mass density in the"
                    " solidification/melting area is assumed.\n"
                    " Please check your settings.\n"
                    " rho0= %5.3e and rho= %5.3e in zone %s\n",
                    __func__, rho0, rho, z->name);

        cs_real_t  cp = cs_property_get_cell_value(z->elt_ids[0],
                                                   time_step->t_cur,
                                                   solid->cp);

        if (fabs(cp - cp0) > FLT_MIN)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: A uniform value of the Cp property in the"
                    " solidification/melting area is assumed.\n"
                    " Please check your settings.\n"
                    " cp0= %5.3e and cp= %5.3e in zone %s\n",
                    __func__, cp0, cp, z->name);

      }

    } /* solidification/melting zone */

  } /* Loop on volume zones */

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    {
      cs_solidification_binary_alloy_t  *alloy
        = (cs_solidification_binary_alloy_t *)solid->model_context;

      if (solid->options & CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM) {

        if (cs_equation_get_space_scheme(alloy->solute_equation) !=
            CS_SPACE_SCHEME_CDOFB)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid space scheme for equation %s\n",
                    __func__, cs_equation_get_name(alloy->solute_equation));

        cs_equation_add_build_hook(alloy->solute_equation,
                                   NULL,                    /* hook context */
                                   _fb_solute_source_term); /* hook function */

        /* Store the pointer to the current face temperature values */

        alloy->temp_faces =
          cs_equation_get_face_values(solid->thermal_sys->thermal_eq, false);

      } /* CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM */

      /* One assumes that all the alloy mixture is liquid thus C_l = C_bulk */

      const size_t  csize = sizeof(cs_real_t)*quant->n_cells;

      memcpy(alloy->c_l_cells, alloy->c_bulk->val, csize);

      /* Set the previous iterate before calling update functions */

      memcpy(alloy->tk_bulk, solid->temperature->val, csize);
      memcpy(alloy->ck_bulk, alloy->c_bulk->val, csize);

      if (alloy->c_l_faces != NULL) {
        cs_real_t  *c_bulk_faces =
          cs_equation_get_face_values(alloy->solute_equation, false);
        memcpy(alloy->c_l_faces, c_bulk_faces,quant->n_faces*sizeof(cs_real_t));
      }

      if (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY)
        _compute_enthalpy(quant,
                          time_step->t_cur,        /* t_eval */
                          solid->temperature->val, /* temperature */
                          solid->g_l_field->val,   /* liquid fraction */
                          alloy->t_eut,            /* temp_ref */
                          solid->latent_heat,      /* latent heat coeff. */
                          solid->mass_density,     /* rho */
                          solid->cp,               /* cp */
                          solid->enthalpy->val);   /* computed enthalpy */

    } /* CS_SOLIDIFICATION_MODEL_BINARY_ALLOY */
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    {
      cs_solidification_voller_t  *v_model
        = (cs_solidification_voller_t *)solid->model_context;

      v_model->update_gl(mesh, connect, quant, time_step);

      if ( (solid->post_flag & CS_SOLIDIFICATION_POST_ENTHALPY) ||
           (solid->model == CS_SOLIDIFICATION_MODEL_VOLLER_NL) )
        _compute_enthalpy(quant,
                          time_step->t_cur,        /* t_eval */
                          solid->temperature->val, /* temperature */
                          solid->g_l_field->val,   /* liquid fraction */
                          v_model->t_solidus,      /* temp_ref */
                          solid->latent_heat,      /* latent heat coeff. */
                          solid->mass_density,     /* rho */
                          solid->cp,               /* cp */
                          solid->enthalpy->val);   /* computed enthalpy */
    }
    break;

  case CS_SOLIDIFICATION_MODEL_STEFAN:
    {
      cs_solidification_stefan_t  *s_model
        = (cs_solidification_stefan_t *)solid->model_context;

      /* Temperature has been initialized.
       * Compute a first guess for the liquid fraction knowing that the liquid
       * fractions is a step function w.r.t. the temperature. Thus, one assumes
       * that the liquid fraction is either 0 or 1 after this step. If the
       * temperature is equal to T_ch, one assumes that g_l = 1
       *
       * Initialize source term and reaction term
       */

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
      for (cs_lnum_t c = 0; c < quant->n_cells; c++) {

        if (solid->temperature->val[c] < s_model->t_change) { /* Solid part */
          solid->g_l_field->val[c] = 0;
          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_SOLID;
        }
        else { /* Liquid part */
          solid->g_l_field->val[c] = 1;
          solid->cell_state[c] = CS_SOLIDIFICATION_STATE_LIQUID;
        }

        /* No source term and reaction term at the begining. One assumes that
           there is no phase change at the first sub-iteration. */

        solid->thermal_reaction_coef_array[c] = 0;
        solid->thermal_source_term_array[c] = 0;

      } /* Loop on cells */

      /* Now compute the enthalpy knowing the temperature and the liquid
         fraction */

      _compute_enthalpy(quant,
                        time_step->t_cur,        /* t_eval */
                        solid->temperature->val, /* temperature */
                        solid->g_l_field->val,   /* liquid fraction */
                        s_model->t_change,       /* temp_ref */
                        solid->latent_heat,      /* latent heat coeff. */
                        solid->mass_density,     /* rho */
                        solid->cp,               /* cp */
                        solid->enthalpy->val);   /* computed enthalpy */

    } /* Stefan model */
    break;

  default:
    break; /* Nothing to do */

  } /* Switch on model */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the solidification module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_compute(const cs_mesh_t              *mesh,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_time_step_t         *time_step)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  switch (solid->model) {

  case CS_SOLIDIFICATION_MODEL_BINARY_ALLOY:
    _default_binary_coupling(mesh, connect, quant, time_step);
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87:
    _voller_prakash_87(mesh, connect, quant, time_step);
    break;

  case CS_SOLIDIFICATION_MODEL_VOLLER_NL:
    _voller_non_linearities(mesh, connect, quant, time_step);
    break;

  case CS_SOLIDIFICATION_MODEL_STEFAN:
      _stefan_thermal_non_linearities(mesh, connect, quant, time_step);
      break;

  default:
    break; /* Nothing else to do */

  } /* Switch on model */

  /* Solve the Navier-Stokes system */

  if ((solid->options & CS_SOLIDIFICATION_NO_VELOCITY_FIELD) == 0)
    /* The Navier-Stokes is not solved when the frozen field is set */
    cs_navsto_system_compute(mesh, connect, quant, time_step);

  /* Perform the monitoring */

  if (solid->verbosity > 0)
    _do_monitoring(quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the solidification module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL)
    return;

  /* Estimate the number of values to output */

  int  n_output_values = CS_SOLIDIFICATION_N_STATES - 1;
  if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {
    n_output_values += 1;

    if (solid->post_flag & CS_SOLIDIFICATION_POST_SEGREGATION_INDEX)
      n_output_values += 1;

  }

  if (solid->post_flag & CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE)
    n_output_values += 1;

  /* Compute the output values */

  cs_real_t  *output_values = NULL;
  BFT_MALLOC(output_values, n_output_values, cs_real_t);
  memset(output_values, 0, n_output_values*sizeof(cs_real_t));

  int n_output_states = (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) ?
    CS_SOLIDIFICATION_N_STATES : CS_SOLIDIFICATION_N_STATES - 1;
  for (int i = 0; i < n_output_states; i++)
    output_values[i] = solid->state_ratio[i];

  n_output_values = n_output_states;

  if (solid->post_flag & CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE) {

    const cs_real_t  *gl = solid->g_l_field->val;

    cs_real_t  integr = 0;
    for (cs_lnum_t i = 0; i < quant->n_cells; i++) {
      if (connect->cell_flag[i] & CS_FLAG_SOLID_CELL)
        continue;
      integr += (1 - gl[i])*quant->cell_vol[i];
    }

    /* Parallel reduction */

    cs_parall_sum(1, CS_REAL_TYPE, &integr);

    output_values[n_output_values] = integr/quant->vol_tot;
    n_output_values++;

  }

  if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *alloy
      = (cs_solidification_binary_alloy_t *)solid->model_context;
    assert(alloy != NULL);

    const cs_real_t  *c_bulk = alloy->c_bulk->val;

    if (solid->post_flag & CS_SOLIDIFICATION_POST_SEGREGATION_INDEX) {

      const cs_real_t  inv_cref = 1./alloy->ref_concentration;

      cs_real_t  si = 0;
      for (cs_lnum_t i = 0; i < quant->n_cells; i++) {
        if (connect->cell_flag[i] & CS_FLAG_SOLID_CELL)
          continue;
        double  c = (c_bulk[i] - alloy->ref_concentration)*inv_cref;
        si += c*c*quant->cell_vol[i];
      }

      /* Parallel reduction */

      cs_parall_sum(1, CS_REAL_TYPE, &si);

      output_values[n_output_values] = sqrt(si/quant->vol_tot);
      n_output_values++;

    }

    if (solid->post_flag & CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE) {
      assert(alloy->t_liquidus != NULL);

      /* Compute the value to be sure that it corresponds to the current
         state */

      for (cs_lnum_t i = 0; i < quant->n_cells; i++) {
        if (connect->cell_flag[i] & CS_FLAG_SOLID_CELL)
          alloy->t_liquidus[i] = -999.99; /* no physical meaning */
        else
          alloy->t_liquidus[i] = _get_t_liquidus(alloy, alloy->c_bulk->val[i]);
      }

    }

    if ((solid->post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS) > 0) {

      assert(alloy->t_liquidus != NULL &&
             alloy->cliq_minus_cbulk != NULL &&
             alloy->tbulk_minus_tliq != NULL);

      const cs_real_t  *c_l = alloy->c_l_cells;
      const cs_real_t  *t_bulk = solid->temperature->val;

      /* Compute Cbulk - Cliq */

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        if (connect->cell_flag[c_id] & CS_FLAG_SOLID_CELL)
          continue; /* = 0 by default */

        const cs_real_t  conc = c_bulk[c_id];
        const cs_real_t  temp = t_bulk[c_id];

        alloy->cliq_minus_cbulk[c_id] = c_l[c_id] - conc;
        alloy->tbulk_minus_tliq[c_id] = temp - alloy->t_liquidus[c_id];

      } /* Loop on cells */

    } /* Advanced analysis */

  } /* Binary alloy modelling */

  if (cs_glob_rank_id < 1 && solid->plot_state != NULL)
    cs_time_plot_vals_write(solid->plot_state,
                            ts->nt_cur,
                            ts->t_cur,
                            n_output_values,
                            output_values);

  BFT_FREE(output_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the solidification module.
 *         Prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_post(void                      *input,
                             int                        mesh_id,
                             int                        cat_id,
                             int                        ent_flag[5],
                             cs_lnum_t                  n_cells,
                             cs_lnum_t                  n_i_faces,
                             cs_lnum_t                  n_b_faces,
                             const cs_lnum_t            cell_ids[],
                             const cs_lnum_t            i_face_ids[],
                             const cs_lnum_t            b_face_ids[],
                             const cs_time_step_t      *time_step)
{
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == NULL)
    return;

  cs_solidification_t  *solid = (cs_solidification_t *)input;

  if (cat_id == CS_POST_MESH_PROBES) {

    cs_field_t  *fld = cs_field_by_name_try("liquid_fraction");
    assert(fld != NULL);

    cs_post_write_probe_values(mesh_id,
                               CS_POST_WRITER_ALL_ASSOCIATED,
                               "liquid_fraction",
                               fld->dim,
                               CS_POST_TYPE_cs_real_t,
                               CS_MESH_LOCATION_CELLS,
                               NULL,
                               NULL,
                               fld->val,
                               time_step);

    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

      cs_solidification_binary_alloy_t  *alloy
        = (cs_solidification_binary_alloy_t *)solid->model_context;

      cs_post_write_probe_values(mesh_id,
                                 CS_POST_WRITER_ALL_ASSOCIATED,
                                 "C_l",
                                 1,
                                 CS_POST_TYPE_cs_real_t,
                                 CS_MESH_LOCATION_CELLS,
                                 NULL,
                                 NULL,
                                 alloy->c_l_cells,
                                 time_step);

      if (solid->post_flag & CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE) {
        assert(alloy->t_liquidus != NULL);
        cs_post_write_probe_values(mesh_id,
                                   CS_POST_WRITER_ALL_ASSOCIATED,
                                   "Tliquidus",
                                   1,
                                   CS_POST_TYPE_cs_real_t,
                                   CS_MESH_LOCATION_CELLS,
                                   NULL,
                                   NULL,
                                   alloy->t_liquidus,
                                   time_step);
      }

      if (solid->post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS) {

        cs_post_write_probe_values(mesh_id,
                                   CS_POST_WRITER_ALL_ASSOCIATED,
                                   "delta_cliq_minus_cbulk",
                                   1,
                                   CS_POST_TYPE_cs_real_t,
                                   CS_MESH_LOCATION_CELLS,
                                   NULL,
                                   NULL,
                                   alloy->cliq_minus_cbulk,
                                   time_step);

        cs_post_write_probe_values(mesh_id,
                                   CS_POST_WRITER_ALL_ASSOCIATED,
                                   "delta_tbulk_minus_tliq",
                                   1,
                                   CS_POST_TYPE_cs_real_t,
                                   CS_MESH_LOCATION_CELLS,
                                   NULL,
                                   NULL,
                                   alloy->tbulk_minus_tliq,
                                   time_step);

        if (alloy->eta_coef_array != NULL)
          cs_post_write_probe_values(mesh_id,
                                     CS_POST_WRITER_ALL_ASSOCIATED,
                                     "Cbulk_advection_scaling",
                                     1,
                                     CS_POST_TYPE_cs_real_t,
                                     CS_MESH_LOCATION_CELLS,
                                     NULL,
                                     NULL,
                                     alloy->eta_coef_array,
                                     time_step);

      } /* Advanced analysis */

    } /* Binary alloy model */

  } /* Probes */

  if ((cat_id == CS_POST_MESH_VOLUME) &&
      (ent_flag[0] == 1)) {     /* ent_flag == 1 --> on cells */

    if (solid->cell_state != NULL &&
        (solid->post_flag & CS_SOLIDIFICATION_POST_CELL_STATE)) {

      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_DEFAULT,
                        "cell_state",
                        1,
                        false,  /* interlace */
                        true,   /* true = original mesh */
                        CS_POST_TYPE_int,
                        solid->cell_state, NULL, NULL,
                        time_step);

    }

    if (solid->model == CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

      cs_solidification_binary_alloy_t  *alloy
        = (cs_solidification_binary_alloy_t *)solid->model_context;

      cs_real_t  *wb = cs_cdo_toolbox_get_tmpbuf();

      if (solid->post_flag & CS_SOLIDIFICATION_ADVANCED_ANALYSIS) {

        if (alloy->cliq_minus_cbulk != NULL)
          cs_post_write_var(CS_POST_MESH_VOLUME,
                            CS_POST_WRITER_DEFAULT,
                            "delta_cliq_minus_cbulk",
                            1,
                            false,  /* interlace */
                            true,   /* true = original mesh */
                            CS_POST_TYPE_cs_real_t,
                            alloy->cliq_minus_cbulk, NULL, NULL,
                            time_step);

        if (alloy->tbulk_minus_tliq != NULL)
          cs_post_write_var(CS_POST_MESH_VOLUME,
                            CS_POST_WRITER_DEFAULT,
                            "delta_tbulk_minus_tliq",
                            1,
                            false,  /* interlace */
                            true,   /* true = original mesh */
                            CS_POST_TYPE_cs_real_t,
                            alloy->tbulk_minus_tliq, NULL, NULL,
                            time_step);

        if (alloy->eta_coef_array != NULL)
          cs_post_write_var(CS_POST_MESH_VOLUME,
                            CS_POST_WRITER_DEFAULT,
                            "Cbulk_advection_scaling",
                            1,
                            false,  /* interlace */
                            true,   /* true = original mesh */
                            CS_POST_TYPE_cs_real_t,
                            alloy->eta_coef_array, NULL, NULL,
                            time_step);

      } /* Advanced analysis */

      if (solid->post_flag & CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE) {

        if (alloy->t_liquidus != NULL)
          cs_post_write_var(CS_POST_MESH_VOLUME,
                            CS_POST_WRITER_DEFAULT,
                            "T_liquidus",
                            1,
                            false,  /* interlace */
                            true,   /* true = original mesh */
                            CS_POST_TYPE_cs_real_t,
                            alloy->t_liquidus, NULL, NULL,
                            time_step);

      }

      if (solid->post_flag & CS_SOLIDIFICATION_POST_CBULK_ADIM) {

        const cs_real_t  inv_cref = 1./alloy->ref_concentration;
        const cs_real_t  *c_bulk = alloy->c_bulk->val;

        for (cs_lnum_t i = 0; i < n_cells; i++)
          wb[i] = (c_bulk[i] - alloy->ref_concentration)*inv_cref;

        cs_post_write_var(CS_POST_MESH_VOLUME,
                          CS_POST_WRITER_DEFAULT,
                          "C_bulk_adim",
                          1,
                          false,  /* interlace */
                          true,   /* true = original mesh */
                          CS_POST_TYPE_cs_real_t,
                          wb, NULL, NULL,
                          time_step);

      } /* CS_SOLIDIFICATION_POST_CBULK_ADIM */

      if (solid->post_flag & CS_SOLIDIFICATION_POST_CLIQ)
        cs_post_write_var(CS_POST_MESH_VOLUME,
                          CS_POST_WRITER_DEFAULT,
                          "C_l",
                          1,
                          false,  /* interlace */
                          true,   /* true = original mesh */
                          CS_POST_TYPE_cs_real_t,
                          alloy->c_l_cells, NULL, NULL,
                          time_step);

    } /* Binary alloy model */

  } /* volume_mesh + on cells */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
