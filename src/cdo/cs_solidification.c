/*============================================================================
 * Handle the solidification module with CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_navsto_system.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solidification.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_SOLIDIFICATION_DBG     0


#define CS_STATE_SOLID   0
#define CS_STATE_LIQUID  1
#define CS_STATE_MUSHY   2
#define CS_N_STATES      3

/*============================================================================
 * Static variables
 *============================================================================*/

static cs_real_t  cs_solidification_forcing_eps  = 1e-3;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer associated to a solidification model aiming at
 *         updating/initializing the solidification variables/properties
 *         dedicated to the model
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_solidification_update_t)(const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_time_step_t        *ts,
                             bool                         cur2prev);


/* Structure storing physical parameters related to a choice of solidification
   modelling */


/* Voller and Prakash model "A fixed grid numerical modelling methodology for
 * convection-diffusion mushy region phase-change problems" Int. J. Heat
 * Transfer, 30 (8), 1987.
 * No tracer. Only physical constants describing the solidification process are
 * used.
 */

typedef struct {

  /* Physical parameters to specify the law of variation of the liquid fraction
   * with respect to the temperature
   *
   * gl(T) = 1 if T > t_liquidus and gl(T) = 0 if T < t_solidus
   * Otherwise:
   * gl(T) = (T - t_solidus)/(t_liquidus - t_solidus)
   */

  cs_real_t        t_solidus;
  cs_real_t        t_liquidus;

  /* Physical parameter for computing the source term in the energy equation
   * Latent heat between the liquid and solid phase
   */
  cs_real_t        latent_heat;

  /* Porous media like reaction term in the momentum equation:
   * F(u) = forcing_coef * (1- gl)^2/(gl^3 + forcing_eps) * u
   */
  cs_real_t        forcing_coef;

} cs_solidification_voller_t;

typedef struct {

  /* Parameters for the Boussinesq approximation in the momentum equation
   * related to the variation of concentration alloy
   */
  cs_real_t        dilatation_coef;
  cs_real_t        ref_concentration;

  /* Physical parameter for computing the source term in the energy equation
   * Latent heat between the liquid and solid phase
   */
  cs_real_t        latent_heat;

  /* Porous media like reaction term in the momentum equation:
   * F(u) = forcing_coef * (1- gl)^2/(gl^3 + forcing_eps) * u
   */
  cs_real_t        forcing_coef;

  /* Phase diagram features */
  /* ---------------------- */

  /* Liquidus slope \frac{\partial g_l}{\partial C} */
  cs_real_t        m_l;

  /* Physical parameters */
  cs_real_t        kp;          /* distribution coefficient */

  /* Alloy features */
  /* -------------- */

  /* The variable related to this equation in the concentration mixture C
   * C = gs*Cs + gl*Cl where gs + gl = 1
   * gl is the liquid fraction and gs the solid fraction
   * Cs = kp * Cl
   *
   * --> C = (gs*kp + gl)*Cl
   */
  cs_equation_t   *tracer_equation;

  /* Tracer concentration in the liquid phase (only if */
  cs_field_t      *c_l_field;

} cs_solidification_binary_alloy_t;


/* Structure storing the set of parameters/structures related to the
 * solidification module */

struct _solidification_t {

  cs_flag_t        model;       /* Modelling for the solidifcation module */
  cs_flag_t        options;     /* Flag dedicated to general options to handle
                                 * the solidification module*/
  cs_flag_t        post_flag;   /* Flag dedicated to the post-processing
                                 * of the solidifcation module */

  /* Mass density of the liquid/solid media */
  cs_property_t   *mass_density;

  /* Advection field (velocity field arising from the Navier-Stokes system) */
  cs_adv_field_t  *adv_field;

  /* Liquid fraction of the mixture */
  /* ------------------------------ */

  cs_field_t      *gl_field;   /* field storing the values of the liquid at
                                  each cell */
  cs_property_t   *gl;         /* liquid fraction property */

  /* array storing the state (solid, mushy, liquid) for each cell */
  int             *cell_state;

  /* Monitoring related to this module */
  cs_real_t        state_ratio[CS_N_STATES];
  cs_gnum_t        n_g_cells[CS_N_STATES];

  /* Quantities related to the energy equation */
  /* ----------------------------------------- */

  /* Fields associated to this module */
  cs_field_t      *temperature;

  /* A reaction term and source term are introduced in the thermal model */
  cs_property_t   *thermal_reaction_coef;
  cs_real_t       *thermal_reaction_coef_array;
  cs_real_t       *thermal_source_term_array;

  /* Additional settings related to the choice of solidification modelling */
  void             *model_context;

  /* A reaction term is introduced in the momentum equation. This terms tends to
   * a huge number when the liquid fraction tends to 0 in order to penalize
   * the velocity to zero when the whole cell is solid
   */
  cs_real_t       *forcing_mom_array; /* values of the forcing reaction
                                         coefficient in each cell */
  cs_property_t   *forcing_mom;

  /* Function pointer related to the way of updating the model */
  cs_solidification_update_t    *update;

};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_module[] =
  " Stop execution.\n"
  " The structure related to the solidifcation module is empty.\n"
  " Please check your settings.\n";

static cs_solidification_t  *cs_solidification_structure = NULL;

/*============================================================================
 * Private static inline function prototypes
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the structure dedicated to the management of the
 *         solidifcation module
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
  solid->model = 0;
  solid->options = 0;
  solid->post_flag = 0;

  /* Property */
  solid->mass_density = NULL;

  /* Advection field arising from the (Navier-)Stokes model */
  solid->adv_field = NULL;

  /* Quantities related to the liquid fraction */
  solid->gl = NULL;
  solid->gl_field = NULL;

  /* State related to each cell */
  solid->cell_state = NULL;

  /* Monitoring */
  for (int i = 0; i < CS_N_STATES; i++) solid->n_g_cells[i] = 0;
  for (int i = 0; i < CS_N_STATES; i++) solid->state_ratio[i] = 0;

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

  /* Function pointer to update the model */
  solid->update = NULL;

  return solid;
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
  cs_adjacency_t  *c2f = connect->c2f;
  cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  cs_solidification_t  *solid = cs_solidification_structure;
  cs_real_t  *face_velocity = cs_xdef_get_array(solid->adv_field->definition);

  /* List of solid cells */
  cs_lnum_t  *solid_cells = NULL;
  BFT_MALLOC(solid_cells, solid->n_g_cells[CS_STATE_SOLID], cs_lnum_t);

  cs_lnum_t  ii = 0;
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
    if (solid->cell_state[c_id] == CS_STATE_SOLID) {
      solid_cells[ii++] = c_id;

      /* Kill the advection field for each face attached to a solid cell */
      for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {
        cs_real_t  *_vel_f = face_velocity + 3*c2f->ids[j];
        _vel_f[0] = _vel_f[1] = _vel_f[2] = 0.;
        }

    }
  } /* Loop on cells */

  assert((cs_gnum_t)ii == solid->n_g_cells[CS_STATE_SOLID]);
  cs_real_t  zero_velocity[3] = {0, 0, 0};
  cs_equation_enforce_by_cell_selection(mom_eqp,
                                        solid->n_g_cells[CS_STATE_SOLID],
                                        solid_cells,
                                        zero_velocity,
                                        NULL);

  BFT_FREE(solid_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the liquid fraction and its related quantities
 *         This corresponds to the Voller and Prakash (87)
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

static void
_update_liquid_fraction_voller(const cs_mesh_t             *mesh,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_time_step_t        *ts,
                               bool                         cur2prev)
{
  CS_UNUSED(mesh);

  cs_solidification_t  *solid = cs_solidification_structure;
  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;

  /* Sanity checks */
  assert(solid->temperature != NULL);
  assert(v_model != NULL);

  if (cur2prev)
    cs_field_current_to_previous(solid->gl_field);

  cs_real_t  *g_l = solid->gl_field->val;
  cs_real_t  *temp = solid->temperature->val;
  assert(temp != NULL);

  /* 1./(t_liquidus - t_solidus) = \partial gl/\partial Temp */
  const cs_real_t  dgl_ov_dtemp = 1./(v_model->t_liquidus - v_model->t_solidus);
  const cs_real_t  inv_forcing_eps = 1./cs_solidification_forcing_eps;

  for (int i = 0; i < CS_N_STATES; i++) solid->n_g_cells[i] = 0;

  /* Retrieve the value of the mass density (assume to be uniform) */
  assert(cs_property_is_uniform(solid->mass_density));
  const cs_real_t  cell_rho = cs_property_get_cell_value(0, ts->t_cur,
                                                         solid->mass_density);
  const cs_real_t  gl_steep =
    cell_rho*v_model->latent_heat*dgl_ov_dtemp/ts->dt[0];

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Update the liquid fraction
     * Update the source term and the reaction coefficient for the thermal
     * system which are arrays */
    if (temp[c_id] < v_model->t_solidus) {

      g_l[c_id] = 0;
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;

      solid->cell_state[c_id] = CS_STATE_SOLID;
      solid->n_g_cells[CS_STATE_SOLID] += 1;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */
      solid->forcing_mom_array[c_id] = v_model->forcing_coef*inv_forcing_eps;

    }
    else if (temp[c_id] > v_model->t_liquidus) {

      g_l[c_id] = 1;
      solid->thermal_reaction_coef_array[c_id] = 0;
      solid->thermal_source_term_array[c_id] = 0;

      solid->n_g_cells[CS_STATE_LIQUID] += 1;
      solid->cell_state[c_id] = CS_STATE_LIQUID;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */
      solid->forcing_mom_array[c_id] = 0;

    }
    else { /* Mushy zone */

      const cs_real_t  glc = (temp[c_id] - v_model->t_solidus) * dgl_ov_dtemp;

      g_l[c_id] = glc;
      solid->thermal_reaction_coef_array[c_id] = gl_steep;
      solid->thermal_source_term_array[c_id] =
        gl_steep*temp[c_id]*quant->cell_vol[c_id];

      solid->cell_state[c_id] = CS_STATE_MUSHY;
      solid->n_g_cells[CS_STATE_MUSHY] += 1;

      /* Update the forcing coefficient treated as a property for a reaction
         term in the momentum eq. */
      const cs_real_t  glm1 = 1 - glc;
      solid->forcing_mom_array[c_id] = v_model->forcing_coef *
        glm1*glm1/(glc*glc*glc + cs_solidification_forcing_eps);

    }

  } /* Loop on cells */

  /* At this stage, the number of solid cells is a local count
   * Set the enforcement of the velocity for solid cells */
  if (solid->n_g_cells[CS_STATE_SOLID] > 0)
    _enforce_solid_cells(connect, quant);

  /* Parallel synchronization of the number of cells in each state */
  cs_parall_sum(CS_N_STATES, CS_GNUM_TYPE, solid->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the monitoring dedicated to the solidifcation module
 *
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_do_monitoring(const cs_cdo_quantities_t   *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  assert(solid->temperature != NULL);

  for (int i = 0; i < CS_N_STATES; i++) solid->state_ratio[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  vol_c = quant->cell_vol[c_id];

    if (solid->cell_state[c_id] == CS_STATE_SOLID)
      solid->state_ratio[CS_STATE_SOLID] += vol_c;
    else if (solid->cell_state[c_id] == CS_STATE_LIQUID)
      solid->state_ratio[CS_STATE_LIQUID] += vol_c;
    else if (solid->cell_state[c_id] == CS_STATE_MUSHY)
      solid->state_ratio[CS_STATE_MUSHY] += vol_c;

  } /* Loop on cells */

  /* Finalize the monitoring step*/
  cs_parall_sum(CS_N_STATES, CS_REAL_TYPE, solid->state_ratio);
  const double  inv_voltot = 1./quant->vol_tot;
  for (int i = 0; i < CS_N_STATES; i++) solid->state_ratio[i] *= inv_voltot;

  cs_log_printf(CS_LOG_DEFAULT,
                "### Solidification monitoring: liquid/mushy/solid states\n");
  cs_log_printf(CS_LOG_DEFAULT, "  * Liquid: %5.2f\%% for %9lu cells;\n",
                100*solid->state_ratio[CS_STATE_LIQUID],
                solid->n_g_cells[CS_STATE_LIQUID]);
  cs_log_printf(CS_LOG_DEFAULT, "  * Mushy:  %5.2f\%% for %9lu cells; \n",
                100*solid->state_ratio[CS_STATE_MUSHY],
                solid->n_g_cells[CS_STATE_MUSHY]);
  cs_log_printf(CS_LOG_DEFAULT, "  * Solid:  %5.2f\%% for %9lu cells; \n",
                100*solid->state_ratio[CS_STATE_SOLID],
                solid->n_g_cells[CS_STATE_SOLID]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term for the momentum equation arising from the
 *         Boussinesq approximation. This relies on the prototype associated
 *         to the generic function pointer \ref cs_dof_function_t
 *
 * \param[in]      n_elts   number of elements to consider
 * \param[in]      elt_ids  list of elements ids
 * \param[in]      compact  true:no indirection, false:indirection for retval
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] retval   result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_solidification_boussinesq_source_term(cs_lnum_t            n_elts,
                                       const cs_lnum_t     *elt_ids,
                                       bool                 compact,
                                       void                *input,
                                       cs_real_t           *retval)
{
  /* Sanity checks */
  assert(retval != NULL && input != NULL);

  /* input is a pointer to a structure */
  const cs_source_term_boussinesq_t  *bq = (cs_source_term_boussinesq_t *)input;

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    cs_lnum_t  id = (elt_ids == NULL) ? i : elt_ids[i]; /* cell_id */
    cs_lnum_t  r_id = compact ? i : id;                 /* position in retval */
    cs_real_t  *_r = retval + 3*r_id;

    /* Thermal effect */
    cs_real_t  bq_coef = -bq->beta*(bq->var[id]-bq->var0);

    for (int k = 0; k < 3; k++)
      _r[k] = bq->rho0 * bq_coef * bq->g[k];

  } /* Loop on elements */
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
 * \brief  Activate the solidification module
 *
 * \param[in]  model            type of modelling
 * \param[in]  options          flag to handle optional parameters
 * \param[in]  post_flag        predefined post-processings
 * \param[in]  boundaries       pointer to the domain boundaries
 * \param[in]  algo_coupling    algorithm used for solving the NavSto system
 * \param[in]  ns_option        option flag for the Navier-Stokes system
 * \param[in]  ns_post_flag     predefined post-processings for Navier-Stokes
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_activate(cs_solidification_model_t      model,
                           cs_flag_t                      options,
                           cs_flag_t                      post_flag,
                           const cs_boundary_t           *boundaries,
                           cs_navsto_param_coupling_t     algo_coupling,
                           cs_flag_t                      ns_option,
                           cs_flag_t                      ns_post_flag)
{
  if (model < 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid modelling. Model = %d\n", __func__, model);

  /* Allocate an empty structure */
  cs_solidification_t  *solid = _solidification_create();

  /* Set members of the structure according to the given settings */
  solid->model = model;
  solid->options = options;
  solid->post_flag = post_flag;

  /* Activate and default settings for the Navier-Stokes module */
  /* ---------------------------------------------------------- */

  cs_navsto_param_model_t  ns_model = CS_NAVSTO_MODEL_SOLIDIFICATION_BOUSSINESQ;
  if (model & CS_SOLIDIFICATION_MODEL_STOKES)
    ns_model |= CS_NAVSTO_MODEL_STOKES;
  else if (model & CS_SOLIDIFICATION_MODEL_NAVIER_STOKES)
    ns_model |= CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES;

  /* Activate the Navier-Stokes module */
  cs_navsto_system_t  *ns = cs_navsto_system_activate(boundaries,
                                                      ns_model,
                                                      algo_coupling,
                                                      ns_option,
                                                      ns_post_flag);

  solid->mass_density = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  assert(solid->mass_density != NULL);

  solid->adv_field = ns->adv_field;

  /* Activate and default settings for the thermal module */
  /* ---------------------------------------------------- */

  cs_flag_t  thm_num = 0, thm_post = 0;
  cs_flag_t  thm_model = CS_THERMAL_MODEL_NAVSTO_VELOCITY;

  if (model & CS_SOLIDIFICATION_MODEL_USE_TEMPERATURE)
    thm_model |= CS_THERMAL_MODEL_USE_TEMPERATURE;
  else if (model & CS_SOLIDIFICATION_MODEL_USE_ENTHALPY)
    thm_model |= CS_THERMAL_MODEL_USE_ENTHALPY;
  else
    thm_model |= CS_THERMAL_MODEL_USE_TEMPERATURE; /* by default */

  cs_thermal_system_activate(thm_model, thm_num, thm_post);

  if (thm_model & CS_THERMAL_MODEL_USE_TEMPERATURE) {

    /* Add reaction property for the temperature equation */
    solid->thermal_reaction_coef = cs_property_add("thermal_reaction_coef",
                                                   CS_PROPERTY_ISO);

    cs_equation_param_t  *th_eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);

    cs_equation_add_reaction(th_eqp, solid->thermal_reaction_coef);

  }

  /* Add properties related to this module */
  solid->forcing_mom = cs_property_add("forcing_momentum_coef",
                                       CS_PROPERTY_ISO);
  solid->gl = cs_property_add("liquid_fraction", CS_PROPERTY_ISO);

  /* Allocate the structure storing the modelling context/settings */

  if (solid->model & CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87) {

    cs_solidification_voller_t  *v_model = NULL;
    BFT_MALLOC(v_model, 1, cs_solidification_voller_t);
    solid->model_context = (void *)v_model;
    solid->update = _update_liquid_fraction_voller;

  }

  else if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *b_model = NULL;
    BFT_MALLOC(b_model, 1, cs_solidification_binary_alloy_t);
    solid->model_context = (void *)b_model;
    solid->update = NULL;       /* TODO */

  }

  /* Set the global pointer */
  cs_solidification_structure = solid;

  return solid;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value of the epsilon parameter used in the forcing term
 *         of the momemtum equation
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
 * \brief  Set the main physical parameters which described the Voller and
 *         Prakash modelling
 *
 * \param[in]  t_solidus      solidus temperature (in K)
 * \param[in]  t_liquidus     liquidus temperatur (in K)
 * \param[in]  latent_heat    latent heat
 * \param[in]  forcing_coef   (< 0) coefficient in the reaction term to reduce
 *                            the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_voller_model(cs_real_t    t_solidus,
                                   cs_real_t    t_liquidus,
                                   cs_real_t    latent_heat,
                                   cs_real_t    forcing_coef)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  if ((solid->model & CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87) == 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Voller and Prakash model not declared during the"
              " activation of the solidification module.\n"
              " Please check your settings.", __func__);

  cs_solidification_voller_t  *v_model
    = (cs_solidification_voller_t *)solid->model_context;
  assert(v_model != NULL);

  /* Model parameters */
  v_model->t_solidus = t_solidus;
  v_model->t_liquidus = t_liquidus;
  v_model->latent_heat = latent_heat;
  v_model->forcing_coef = forcing_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the main physical parameters which described a solidification
 *         process with a binary alloy
 *         Add a tracer equation to simulation the conv/diffusion of the alloy
 *         ratio between the two components of the alloy
 *
 * \param[in]  name          name of the binary alloy
 * \param[in]  varname       name of the unknown related to the tracer eq.
 * \param[in]  conc0         reference concentration
 * \param[in]  beta          value of the dilatation coefficient w.r.t. concen
 * \param[in]  kp            value of the distribution coefficient
 * \param[in]  m_l           liquidus slope for the tracer
 * \param[in]  latent_heat   latent heat
 * \param[in]  forcing_coef  (< 0) coefficient in the reaction term to reduce
 *                           the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_binary_alloy_model(const char     *name,
                                         const char     *varname,
                                         cs_real_t       conc0,
                                         cs_real_t       beta,
                                         cs_real_t       kp,
                                         cs_real_t       m_l,
                                         cs_real_t       latent_heat,
                                         cs_real_t       forcing_coef)
{
  cs_solidification_t  *solid = cs_solidification_structure;
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  cs_solidification_binary_alloy_t  *b_model
    = (cs_solidification_binary_alloy_t *)solid->model_context;

  /* Sanity checks */
  assert(name != NULL && varname != NULL);
  assert(solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY);
  assert(b_model != NULL);

  b_model->tracer_equation = cs_equation_add(name, varname,
                                             CS_EQUATION_TYPE_SOLIDIFICATION,
                                             1,
                                             CS_PARAM_BC_HMG_NEUMANN);

  b_model->c_l_field = NULL;

  /* Set the main physical parameters */
  b_model->dilatation_coef = beta;
  b_model->ref_concentration = conc0;
  b_model->kp = kp;
  b_model->m_l = m_l;
  b_model->latent_heat = latent_heat;
  b_model->forcing_coef = forcing_coef;
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

  if (solid->model & CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87) {
    cs_solidification_voller_t  *v_model
      = (cs_solidification_voller_t *)solid->model_context;

    BFT_FREE(v_model);

  } /* Voller and Prakash modelling */

  if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {
    cs_solidification_binary_alloy_t  *b_model
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    BFT_FREE(b_model);

  } /* Binary alloy modelling */

  BFT_FREE(solid->thermal_reaction_coef_array);
  BFT_FREE(solid->thermal_source_term_array);
  BFT_FREE(solid->forcing_mom_array);

  BFT_FREE(solid->cell_state);

  BFT_FREE(solid);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the Solidification module
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_setup(void)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  const int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_CDO;
  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const int  c_loc_id = cs_mesh_location_get_id_by_name("cells");

  /* Add a field for the liquid fraction */
  solid->gl_field = cs_field_create("liquid_fraction",
                                    field_mask,
                                    c_loc_id,
                                    1,
                                    true); /* has_previous */

  cs_field_set_key_int(solid->gl_field, log_key, 1);
  cs_field_set_key_int(solid->gl_field, post_key, 1);

  /* Add a reaction term to the momentum equation */
  cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  assert(mom_eqp != NULL);

  cs_equation_add_reaction(mom_eqp, solid->forcing_mom);

  /* Add default post-processing related to the solidifcation module */
  cs_post_add_time_mesh_dep_output(cs_solidification_extra_post, solid);

  /* Model-specific part */

  if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *b_model
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    b_model->c_l_field = cs_field_create("alloy_liquid_conc",
                                         field_mask,
                                         c_loc_id,
                                         1,
                                         true); /* has_previous */

    cs_field_set_key_int(b_model->c_l_field, log_key, 1);
    cs_field_set_key_int(b_model->c_l_field, post_key, 1);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the solidification
 *         module
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_finalize_setup(const cs_cdo_connect_t       *connect,
                                 const cs_cdo_quantities_t    *quant)
{
  CS_UNUSED(connect);

  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  /* Retrieve the field associated to the temperature */
  solid->temperature = cs_field_by_name("temperature");

  /* Define the liquid fraction */
  cs_property_def_by_field(solid->gl, solid->gl_field);

  /* Initially one assumes that all is liquid */
  cs_field_set_values(solid->gl_field, 1.);

  BFT_MALLOC(solid->cell_state, quant->n_cells, int);
# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_cells; i++)
    solid->cell_state[i] = CS_STATE_LIQUID;

  /* Add the boussinesq source term in the momentum equation */
  cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  cs_physical_constants_t  *phy_constants = cs_get_glob_physical_constants();
  cs_property_t  *mass_density = solid->mass_density;
  assert(mass_density != NULL && mom_eq != NULL);

  if (solid->model & CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87) {

    /* Define the metadata to build a Boussinesq source term related to the
     * temperature. This structure is allocated here but the lifecycle is
     * managed by the cs_thermal_system_t structure */
    cs_source_term_boussinesq_t  *thm_bq =
      cs_thermal_system_add_boussinesq_source_term(phy_constants->gravity,
                                                   mass_density->ref_value);

    cs_dof_func_t  *func = _solidification_boussinesq_source_term;
    cs_equation_add_source_term_by_dof_func(mom_eqp,
                                            NULL, /* = all cells */
                                            cs_flag_primal_cell,
                                            func,
                                            thm_bq);

  }
  else if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *b_model
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    /* TODO */

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Case with tracer(s) is not handled yet.", __func__);

  /* Define the forcing term acting as a reaction term in the momentum equation
     This term is related to the liquid fraction */
  BFT_MALLOC(solid->forcing_mom_array, quant->n_cells, cs_real_t);
  memset(solid->forcing_mom_array, 0, quant->n_cells*sizeof(cs_real_t));

  cs_property_def_by_array(solid->forcing_mom,
                           cs_flag_primal_cell,
                           solid->forcing_mom_array,
                           false, /* definition is owner ? */
                           NULL); /* no index */

  /* Define the reaction coefficient and the source term for the temperature
     equation */
  if (solid->thermal_reaction_coef != NULL) {

    size_t  array_size = quant->n_cells*sizeof(cs_real_t);
    BFT_MALLOC(solid->thermal_reaction_coef_array, quant->n_cells, cs_real_t);
    memset(solid->thermal_reaction_coef_array, 0, array_size);

    cs_property_def_by_array(solid->thermal_reaction_coef,
                             cs_flag_primal_cell,
                             solid->thermal_reaction_coef_array,
                             false, /* definition is owner ? */
                             NULL); /* no index */

    BFT_MALLOC(solid->thermal_source_term_array, quant->n_cells, cs_real_t);
    memset(solid->thermal_source_term_array, 0, array_size);

    cs_equation_param_t  *thm_eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);
    cs_equation_add_source_term_by_array(thm_eqp,
                                         NULL,   /* all cells selected */
                                         cs_flag_primal_cell,
                                         solid->thermal_source_term_array,
                                         false,  /* definition is owner ? */
                                         NULL);  /* no index */

  }

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

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the solidification module\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", h1_sep);

  cs_log_printf(CS_LOG_SETUP, "  * Solidification | Model:");
  if (solid->model & CS_SOLIDIFICATION_MODEL_STOKES)
    cs_log_printf(CS_LOG_SETUP, "Stokes");
  else if (solid->model & CS_SOLIDIFICATION_MODEL_NAVIER_STOKES)
    cs_log_printf(CS_LOG_SETUP, "Navier-Stokes");
  cs_log_printf(CS_LOG_SETUP, "\n");

  cs_log_printf(CS_LOG_SETUP, "  * Solidification | Model:");
  if (solid->model & CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87) {

    cs_solidification_voller_t  *v_model
      = (cs_solidification_voller_t *)solid->model_context;

    cs_log_printf(CS_LOG_SETUP, "Voller-Prakash (1987)\n");
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Tliq: %5.3e; Tsol: %5.3e",
                  v_model->t_liquidus, v_model->t_solidus);
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Latent heat: %5.3e\n",
                  v_model->latent_heat);
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Forcing coeff: %5.3e\n",
                  v_model->forcing_coef);

  }
  else if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *b_model
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    cs_log_printf(CS_LOG_SETUP, "Binary alloy\n");
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Alloy: %s",
                  cs_equation_get_name(b_model->tracer_equation));
    cs_log_printf(CS_LOG_SETUP,
                  "  * Solidification | Dilatation coef. concentration: %5.3e\n",
                  b_model->dilatation_coef);
    cs_log_printf(CS_LOG_SETUP,
                  "  * Solidification | Reference concentration: %5.3e\n",
                  b_model->ref_concentration);
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Latent heat: %5.3e\n",
                  b_model->latent_heat);
    cs_log_printf(CS_LOG_SETUP, "  * Solidification | Forcing coeff: %5.3e\n",
                  b_model->forcing_coef);

  }

  cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the solidification module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_compute(const cs_mesh_t              *mesh,
                          const cs_time_step_t         *time_step,
                          const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant)
{
  cs_solidification_t  *solid = cs_solidification_structure;

  /* Sanity checks */
  if (solid == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_module));

  if (solid->model & CS_SOLIDIFICATION_MODEL_BINARY_ALLOY) {

    cs_solidification_binary_alloy_t  *b_model
      = (cs_solidification_binary_alloy_t *)solid->model_context;

    cs_equation_solve(mesh, b_model->tracer_equation);
  }

  /* Add equations to be solved at each time step */
  cs_thermal_system_compute(mesh, time_step, connect, quant);

  /* Solve the Navier-Stokes system */
  cs_navsto_system_compute(mesh, time_step, connect, quant);

  /* Update fields and properties which are related to solved variables */
  solid->update(mesh, connect, quant, time_step,
                true); /* operate current to previous ? */

  /* Perform the monitoring */
  _do_monitoring(quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the solidification module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant)
{
  CS_UNUSED(connect);
  CS_UNUSED(quant);

  cs_solidification_t  *solid = cs_solidification_structure;

  if (solid == NULL)
    return;

  /* TODO */
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
  CS_UNUSED(mesh_id);
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == NULL)
    return;

  cs_solidification_t  *solid = (cs_solidification_t *)input;

  if (solid->cell_state != NULL) {

    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_DEFAULT,
                      "cell_state",
                      1,
                      false,  // interlace
                      true,   // true = original mesh
                      CS_POST_TYPE_int,
                      solid->cell_state,
                      NULL,
                      NULL,
                      time_step);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
