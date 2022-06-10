#ifndef __CS_SOLIDIFICATION_H__
#define __CS_SOLIDIFICATION_H__

/*============================================================================
 * Header to handle the solidification module with CDO schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_equation.h"
#include "cs_iter_algo.h"
#include "cs_navsto_param.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @name Flags specifying automatic post-processing for the solidification
 *        module
 * @{
 *
 * \def CS_SOLIDIFICATION_POST_CELL_STATE
 * \brief State related to each cell between (solid, mushy, liquid or eutectic)
 *
 * \def CS_SOLIDIFICATION_POST_ENTHALPY
 * \brief Enthalpy in each cell. By default, only the temperature is
 *        post-processed
 *
 * \def CS_SOLIDIFICATION_POST_CBULK_ADIM
 * \brief Compute and post-process (C_bulk - C_0)/C_0
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated
 *        C_0 is the reference concentration
 *
 * \def CS_SOLIDIFICATION_POST_CLIQ
 * \brief Post-process Cliq the liquid solute distribution (wt %)
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated.
 *
 * \def CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE
 * \brief Activate the (volumic) post-processing of the liquidus temperature
 *        in each cell
 *
 * \def CS_SOLIDIFICATION_POST_SEGREGATION_INDEX
 * \brief Activate the computation and output in the file solidification.dat
 *        for each time step of the segregation index defined by
 *        sqrt( 1/|Domaine| * \int_{Domain} ((C_bulk - C_0)/C_0)**2 )
 *        Only available if the model \ref CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *        is activated
 *
 * \def CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE
 * \brief Activate the computation and output in the file solidification.dat
 *        for each time step of the integral over the computational domain
 *        of the solid fraction divided by the volume of the domain.
 *
 * \def CS_SOLIDIFICATION_ADVANCED_ANALYSIS
 * \brief Activate a set of post-processing (Advanced usage. Only for the
 * understanding of the solidification process)
*/

#define CS_SOLIDIFICATION_POST_CELL_STATE             (1 << 0) /* =   1 */
#define CS_SOLIDIFICATION_POST_ENTHALPY               (1 << 1) /* =   2 */
#define CS_SOLIDIFICATION_POST_CBULK_ADIM             (1 << 2) /* =   4 */
#define CS_SOLIDIFICATION_POST_CLIQ                   (1 << 3) /* =   8 */
#define CS_SOLIDIFICATION_POST_LIQUIDUS_TEMPERATURE   (1 << 4) /* =  16 */
#define CS_SOLIDIFICATION_POST_SEGREGATION_INDEX      (1 << 5) /* =  32 */
#define CS_SOLIDIFICATION_POST_SOLIDIFICATION_RATE    (1 << 6) /* =  64 */
#define CS_SOLIDIFICATION_ADVANCED_ANALYSIS           (1 << 7) /* = 128 */

/*!
 * @name Flags specifying options specific to the solidification module
 * @{
 *
 * \def CS_SOLIDIFICATION_USE_ENTHALPY_VARIABLE
 * \brief The dynamic system of equations is associated with an energy equation
 *        solved using the enthalpy as variable (not fully available).
 *
 * \def CS_SOLIDIFICATION_NO_VELOCITY_FIELD
 * \brief The system of equations does not involve the Navier-Stokes equations.
 *        No velocity is taken into account.
 *
 * \def CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM
 * \brief The solute equation related to the transport of the bulk concentration
 *        is treated with a source term related to an explicit advection of the
 *        quantity (C - Cl). The default behavior is to add a weighting
 *        coefficient to the (implicit) advection term related to the liquid
 *        fraction
 *        This option is related to the LEGACY strategy.
 *
 * \def CS_SOLIDIFICATION_USE_EXTRAPOLATION
 * \brief Use an extrapolation during the computation of different terms
 *        according to the strategy. This extrapolation of variable at time
 *        step n+1 uses values at n and n-1: 2*u^n - u^{n-1}
 *
 * \def CS_SOLIDIFICATION_WITH_PENALIZED_EUTECTIC
 * \brief Option related to the PATH strategy.
 *        Introduce a reaction term and a source term in order to remain on
 *        the eutectic plateau.
 */

#define CS_SOLIDIFICATION_USE_ENTHALPY_VARIABLE             (1 << 0) /*=    1 */
#define CS_SOLIDIFICATION_NO_VELOCITY_FIELD                 (1 << 1) /*=    2 */
#define CS_SOLIDIFICATION_WITH_SOLUTE_SOURCE_TERM           (1 << 2) /*=    4 */
#define CS_SOLIDIFICATION_USE_EXTRAPOLATION                 (1 << 3) /*=    8 */
#define CS_SOLIDIFICATION_WITH_PENALIZED_EUTECTIC           (1 << 4) /*=   16 */

/* Automatically set by the code if user functions are used
 * The following flags are set when calling \ref cs_solidification_set_functions
 *
 * \def CS_SOLIDIFICATION_BINARY_ALLOY_M_FUNC
 * \brief the update of the forcing term (penalization) in the momentum equation
 *        is defined using a user function
 *
 * \def CS_SOLIDIFICATION_BINARY_ALLOY_C_FUNC
 * \brief the update of the liquid concentration of the binary alloy is defined
 *        using a user function
 *
 * \def CS_SOLIDIFICATION_BINARY_ALLOY_G_FUNC
 * \brief the update of the liquid fraction is defined using a user function
 *
 * \def CS_SOLIDIFICATION_BINARY_ALLOY_T_FUNC
 * \brief the update of the thermal source term is defined using a user function
 */
#define CS_SOLIDIFICATION_BINARY_ALLOY_M_FUNC               (1 << 7) /*=  128 */
#define CS_SOLIDIFICATION_BINARY_ALLOY_C_FUNC               (1 << 8) /*=  256 */
#define CS_SOLIDIFICATION_BINARY_ALLOY_G_FUNC               (1 << 9) /*=  512 */
#define CS_SOLIDIFICATION_BINARY_ALLOY_T_FUNC               (1 <<10) /*= 1024 */

/*!
 * @}
 */

/*=============================================================================
 * Structure and type definitions
 *============================================================================*/

/*! \enum cs_solidification_model_t
 *  \brief Type of physical model used to simulate the solidifcation/fusion
 *         process
 *
 * \var CS_SOLIDIFICATION_MODEL_STEFAN
 *      Phase change model without advection field. The phase change is assumed
 *      to be at a given temperature meaning that the liquid fraction is a step
 *      function w.r.t. the temperature.
 *
 * \var CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87
 *      Modelling introduced in Voller and Prakash entitled: "A fixed grid
 *      numerical modelling methodology for convection-diffusion mushy region
 *      phase-change problems" Int. J. Heat Transfer, 30 (8), 1987.  No
 *      tracer. Only physical constants describing the solidification process
 *      are used.
 *
 * \var CS_SOLIDIFICATION_MODEL_VOLLER_NL
 *      Modelling based on \ref CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87 but
 *      the thermal equation and the update of the liquid fraction is non-linear
 *
 * \var CS_SOLIDIFICATION_MODEL_BINARY_ALLOY
 *      The basis is similar to \ref CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87
 *      A tracer equation is added corresponding to the evolution of the bulk
 *      concentration in alloy. This alloy has two chemical constituents hence
 *      the name "binary alloy".
 */

typedef enum {

  CS_SOLIDIFICATION_MODEL_STEFAN,
  CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87,
  CS_SOLIDIFICATION_MODEL_VOLLER_NL,
  CS_SOLIDIFICATION_MODEL_BINARY_ALLOY,

  CS_SOLIDIFICATION_N_MODELS

} cs_solidification_model_t;

/*! \enum cs_solidification_state_t
 *  \brief Kind of state in which a cell or an entity is
 *
 * \var CS_SOLIDIFICATION_STATE_SOLID
 *      Solid state (the liquid fraction is equal to 0)
 *
 * \var CS_SOLIDIFICATION_STATE_MUSHY
 *      Mushy state meaning that the liquid fraction is strictly below 1 and
 *      strictly above 0.
 *
 * \var CS_SOLIDIFICATION_STATE_LIQUID
 *      Liquid state (the liquid fraction is equal to 1)
 *
 * \var CS_SOLIDIFICATION_STATE_EUTECTIC
 *      Eutectic state. Only possible with a \ref
 *      CS_SOLIDIFICATION_MODEL_BINARY_ALLOY model. A rough transition between
 *      the mushy and solid state occurs when an eutectic transition happens.
 */

typedef enum {

  CS_SOLIDIFICATION_STATE_SOLID    = 0,
  CS_SOLIDIFICATION_STATE_MUSHY    = 1,
  CS_SOLIDIFICATION_STATE_LIQUID   = 2,
  CS_SOLIDIFICATION_STATE_EUTECTIC = 3,

  CS_SOLIDIFICATION_N_STATES       = 4,

} cs_solidification_state_t;

/*! \enum cs_solidification_strategy_t
 *  \brief Kind of strategy to use to model the segregation/solidification
 *         process. This implies a setting of functions to update the liquid
 *         fraction, the thermal source terms, the liquid concentration and its
 *         related quantities.
 */

typedef enum {

  CS_SOLIDIFICATION_STRATEGY_LEGACY,
  CS_SOLIDIFICATION_STRATEGY_TAYLOR,
  CS_SOLIDIFICATION_STRATEGY_PATH,

  CS_SOLIDIFICATION_N_STRATEGIES

} cs_solidification_strategy_t;

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
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_solidification_func_t)(const cs_mesh_t             *mesh,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts);

/* Structure storing physical parameters related to a choice of solidification
   modelling */

/*----------------------------------------------------------------------------
 * Stefan model
 *----------------------------------------------------------------------------
 *
 * Neither advection nor segregation is taken into account.
 * The liquid fraction is a step function w.r.t. the temperature.
 */

typedef struct {

  /* Physical parameters to specify the law of variation of the liquid fraction
   * with respect to the temperature
   *
   * gl(T) = 1 if T > T_change
   * Otherwise:
   * gl(T) = 0
   */

  cs_real_t                      t_change;

  /* Function pointers */
  /* ----------------- */

  /* Function to update the liquid fraction */

  cs_solidification_func_t      *update_gl;

  /* Function to update the source term for the thermal equation */

  cs_solidification_func_t      *update_thm_st;

  /* Numerical parameters */
  /* -------------------- */

  int                            n_iter_max;
  double                         max_delta_h;

} cs_solidification_stefan_t;


/* Voller and Prakash model "A fixed grid numerical modelling methodology for
 * convection-diffusion mushy region phase-change problems" Int. J. Heat
 * Transfer, 30 (8), 1987.
 * No tracer. Only physical constants describing the solidification process are
 * used.
 */

/*----------------------------------------------------------------------------
 * Solidification without segregation (Voller & Prakash'87 model)
 *----------------------------------------------------------------------------*/

typedef struct {

  /* Secondary dendrite arm spacing */

  cs_real_t                      s_das;

  /* Physical parameters to specify the law of variation of the liquid fraction
   * with respect to the temperature
   *
   * gl(T) = 1 if T > t_liquidus and gl(T) = 0 if T < t_solidus
   * Otherwise:
   * gl(T) = (T - t_solidus)/(t_liquidus - t_solidus)
   */

  cs_real_t                      t_solidus;
  cs_real_t                      t_liquidus;

  /* Function pointers */
  /* ----------------- */

  /* Function to update the liquid fraction */

  cs_solidification_func_t      *update_gl;

  /* Function to update the source term for the thermal equation */

  cs_solidification_func_t      *update_thm_st;

  /* Numerical parameters */
  /* -------------------- */

  cs_param_nl_algo_t             nl_algo_type;
  cs_iter_algo_t                *nl_algo;

} cs_solidification_voller_t;


/*----------------------------------------------------------------------------
 * Solidification of a binary alloy with segregation (Voller & Prakash'89 model)
 *----------------------------------------------------------------------------*/

typedef struct {

  /* Alloy features */
  /* -------------- */

  /* Reference mixture concentration (used in the Boussinesq approximation and
   * for normalization */

  cs_real_t    ref_concentration;

  /* Phase diagram features for an alloy with the component A and B */
  /* -------------------------------------------------------------- */

  /* Secondary dendrite arm spacing */

  cs_real_t    s_das;

  /* Physical parameters */

  cs_real_t    kp;       /* distribution coefficient */
  cs_real_t    inv_kp;   /* reciprocal of kp */
  cs_real_t    inv_kpm1; /* 1/(kp - 1) */
  cs_real_t    ml;       /* Liquidus slope \frac{\partial g_l}{\partial C} */
  cs_real_t    inv_ml;   /* reciprocal of ml */

  /* Temperature of phase change for the pure material (conc = 0) */

  cs_real_t    t_melt;

  /* Eutectic point: temperature and concentration */

  cs_real_t    t_eut;
  cs_real_t    t_eut_inf;
  cs_real_t    t_eut_sup;

  cs_real_t    c_eut;
  cs_real_t    cs1;
  cs_real_t    dgldC_eut;

  /* The variable related to this equation in the solute concentration of
   * the mixture: c_bulk (c_s in the solid phase and c_l in the liquid phase)
   * c_bulk = gs*c_s + gl*c_l where gs + gl = 1
   * gl is the liquid fraction and gs the solid fraction
   * c_s = kp * c_l (lever rule is assumed up to now)
   *
   * --> c_bulk = (gs*kp + gl)*c_l
   */

  /* Function to update the liquid fraction */

  cs_solidification_func_t     *update_gl;

  /* Function to update the source term for the thermal equation */

  cs_solidification_func_t     *update_thm_st;

  /* Function to update the velocity forcing in the momentum equation */

  cs_solidification_func_t     *update_velocity_forcing;

  /* Function to update c_l in each cell */

  cs_solidification_func_t     *update_clc;

  /* Drive the convergence of the coupled system (solute transport and thermal
   * equation) with respect to the following criteria (taken from Voller and
   * Swaminathan'91)
   *   max_{c\in C} |Temp^(k+1) - Temp^(k)| < delta_tolerance
   *   max_{c\in C} |Cbulk^(k+1) - Cbulk*^(k)| < delta_tolerance
   *   n_iter < n_iter_max
   *
   * eta_relax: add a relaxation in the update of the eta coefficient
   * conc_liq = eta_coef * conc_bulk
   * eta_relax = 0. --> No relaxation (default choice)
   *
   * gliq_relax: idem but for the liquid fraction
   */

  int                              iter;
  int                              n_iter_max;
  double                           delta_tolerance;
  double                           eta_relax;
  double                           gliq_relax;

  /* During the non-linear iteration process one needs:
   *  temp_{n}         --> stored in field->val_pre
   *  temp_{n+1}^k     --> stored in tk_bulk (in this structure)
   *  temp_{n+1}^{k+1} --> stored in field->val
   *
   * Optionally one may consider an extrapolated temperature and bulk
   * concentration
   *  temp_{n+1}^{extrap} = 2*temp_{n} - temp_{n-1}
   *
   * Same thing for the bulk concentration.
   */

  cs_real_t         *tk_bulk;
  cs_real_t         *ck_bulk;
  cs_real_t         *tx_bulk;
  cs_real_t         *cx_bulk;

  /* Solute concentration in the liquid phase
   * 1) array of the last computed values at cells
   * 2) array of the last computed values at faces (interior and border) */

  cs_real_t         *c_l_cells;
  cs_real_t         *c_l_faces;

  /* Temperature values at faces (this is not owned by the structure) */

  const cs_real_t   *temp_faces;

  /* Equation for the solute transport and related quantities */

  cs_equation_t     *solute_equation;
  cs_field_t        *c_bulk;

  /* Diffusion coefficient for the solute in the liquid phase
   * diff_pty_val = rho * g_l * diff_coef */

  cs_real_t          diff_coef;
  cs_property_t     *diff_pty;
  cs_real_t         *diff_pty_array;

  cs_property_t     *eta_coef_pty;
  cs_real_t         *eta_coef_array;

  /* Optional post-processing arrays */
  /* ------------------------------ */

  /* Liquidus temperature (values at cell centers) */

  cs_real_t         *t_liquidus;

  /* Quantities for advanced analysis */

  cs_real_t         *tbulk_minus_tliq;
  cs_real_t         *cliq_minus_cbulk;

} cs_solidification_binary_alloy_t;

/*----------------------------------------------------------------------------
 * Main structure to manage the solidification process
 *----------------------------------------------------------------------------*/

typedef struct  {

  cs_flag_t        model;       /* Modelling for the solidification module */
  cs_flag_t        options;     /* Flag dedicated to general options to handle
                                 * the solidification module*/
  cs_flag_t        post_flag;   /* Flag dedicated to the post-processing
                                 * of the solidification module */
  int              verbosity;   /* Level of verbosity */

  /* Physical properties common to all models */
  /* ---------------------------------------- */

  /* Mass density of the liquid/solid media */

  cs_property_t   *mass_density;

  /* Reference value for the heat capacity in the solidification/melting area
   * (assumed to be uniform) */

  cs_property_t   *cp;

  /* Viscosity (pointer to the total viscosity which should be equal to the
  *  laminar viscosity since no turbulence modelling is usually taken into
  *  account */

  cs_property_t   *viscosity;

  /* Physical parameter for computing the source term in the energy equation
   * Latent heat between the liquid and solid phase
   */

  cs_real_t        latent_heat;


  /* Liquid fraction of the mixture */
  /* ------------------------------ */

  cs_field_t      *g_l_field;   /* field storing the values of the liquid
                                   fraction at each cell */
  cs_property_t   *g_l;         /* liquid fraction property */

  /* array storing the state (solid, mushy, liquid) for each cell */

  cs_solidification_state_t     *cell_state;

  /* Plot evolution of the solidification process */

  cs_time_plot_t                *plot_state;

  /* Monitoring related to this module */

  cs_real_t        state_ratio[CS_SOLIDIFICATION_N_STATES];
  cs_gnum_t        n_g_cells[CS_SOLIDIFICATION_N_STATES];

  /* Quantities related to the energy equation */
  /* ----------------------------------------- */

  cs_thermal_system_t   *thermal_sys;

  /* Fields associated to this module */

  cs_field_t      *temperature;

  /* Enthalpy (by default this a derived field which can be used to update of
     the liquid fraction) */

  cs_field_t      *enthalpy;


  /* A reaction term and source term are introduced in the thermal model */

  cs_property_t   *thermal_reaction_coef;
  cs_real_t       *thermal_reaction_coef_array;
  cs_real_t       *thermal_source_term_array;

  /* Additional settings related to the choice of solidification modelling */

  void            *model_context;

  /* Strategy to update quantities during the solidification process. These
   * quantities are the liquid fraction, the thermal source term for
   * instance */

  cs_solidification_strategy_t   strategy;

  /* A reaction term is introduced in the momentum equation. This terms tends to
   * a huge number when the liquid fraction tends to 0 in order to penalize
   * the velocity to zero when the whole cell is solid
   */

  cs_real_t       *forcing_mom_array; /* values of the forcing reaction
                                         coefficient in each cell */
  cs_property_t   *forcing_mom;

  /* Porous media like reaction term in the momentum equation:
   *
   * forcing_coef = 180 * visco0 / s_das^2
   * where visco0 is the laminar viscosity and s_das is the secondary
   * dendrite arm spacing
   * F(u) = forcing_coef * (1- gl)^2/(gl^3 + forcing_eps) * u
   */

  cs_real_t        forcing_coef;

  /* First cell associated to a fluid/solid area (i.e. not associated to
   * a permanent solid zone) */

  cs_lnum_t        first_cell;

} cs_solidification_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if solidification module is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_solidification_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main structure to deal with solidification process
 *
 * \return a pointer to a new allocated solidification structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_get_structure(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the level of verbosity for the solidification module
 *
 * \param[in]   verbosity     level of verbosity to set
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_verbosity(int   verbosity);

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
cs_solidification_set_forcing_eps(cs_real_t    forcing_eps);

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
cs_solidification_activate(cs_solidification_model_t      model,
                           cs_flag_t                      options,
                           cs_flag_t                      post_flag,
                           const cs_boundary_t           *boundaries,
                           cs_navsto_param_model_t        ns_model,
                           cs_navsto_param_model_flag_t   ns_model_flag,
                           cs_navsto_param_coupling_t     algo_coupling,
                           cs_navsto_param_post_flag_t    ns_post_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the Stefan model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_stefan_t *
cs_solidification_get_stefan_struct(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the Stefan's model settings
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_stefan_t *
cs_solidification_check_stefan_model(void);

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
                                   cs_real_t    latent_heat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the Voller model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_voller_t *
cs_solidification_get_voller_struct(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the Voller's model settings
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_voller_t *
cs_solidification_check_voller_model(void);

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
                                   cs_real_t    s_das);

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
                                               cs_real_t    latent_heat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the structure defining the binary alloy model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_binary_alloy_t *
cs_solidification_get_binary_alloy_struct(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Sanity checks on the consistency of the settings of the binary alloy
 *         model
 *
 * \return a pointer to the structure
 */
/*----------------------------------------------------------------------------*/

cs_solidification_binary_alloy_t *
cs_solidification_check_binary_alloy_model(void);

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
                                         cs_real_t       s_das);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the strategy to update quantitiess (liquid fraction and
 *         the thermal source term for the two main quantities)
 *
 * \param[in]  strategy     strategy to perform the update of quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_set_strategy(cs_solidification_strategy_t  strategy);

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
                                     cs_solidification_func_t  *thm_st_update);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the solidification module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_solidification_t *
cs_solidification_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the solidification module
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_init_setup(void);

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
                                 const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the solidification module in the log file dedicated to
 *         the setup
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_log_setup(void);

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
                              const cs_time_step_t         *time_step);

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
                          const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the solidification module
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_solidification_extra_op(const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the solidification module.
 *         Prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to an optional structure (here a
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
                             const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLIDIFICATION_H__ */
