#ifndef __CS_GWF_PRIV_H__
#define __CS_GWF_PRIV_H__

/*============================================================================
 * Structures/types related to the groundwater flow module
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

#include "cs_advection_field.h"
#include "cs_equation_system.h"
#include "cs_gwf_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_gwf_darcy_flux_t
 *
 * \brief Structure to handle the Darcy flux
 */

typedef struct {

  /*!
   * \var adv_field
   * Pointer to a \ref cs_adv_field_t structure. Darcy advective flux in the
   * liquid phase. This structure is used to define the advective term in
   * tracer equations.
   *
   * \var flux_location
   * Indicate where the arrays defining the Darcy fluxes are located
   *
   * \var flux_val
   * Array storing the liquid Darcian flux in each location (for instance the
   * dual faces associated to each cell)
   *
   * \var boundary_flux_val
   * Array storing the normal Darcian flux across the boundary of the
   * computational domain for the liquid phase. This is an optional array.
   */

  cs_adv_field_t               *adv_field;
  cs_flag_t                     flux_location;
  cs_real_t                    *flux_val;
  cs_real_t                    *boundary_flux_val;

} cs_gwf_darcy_flux_t;


/*! \struct cs_gwf_saturated_single_phase_t
 *
 * \brief Structure to handle the modelling of a single-phase flows in a porous
 *        media considered as saturated.
 *
 * Several simplifications are operated in this context. Only the liquid phase
 * is taken into account.
 */

typedef struct {

  /*!
   * @name Equation
   * @{
   *
   * \var richards

   * The Richards equation is the governing equation which corresponds to the
   * mass conservation of water. "hydraulic_head" is the associated variable
   * "permeability" is the diffusion property related to this equation.
   */

  cs_equation_t                *richards;

  /*!
   * @}
   * @name Darcy (advection) field
   * @{
   *
   * \var darcy
   * Pointer to a \ref cs_gwf_darcy_flux_t structure. Darcy advective flux in
   * the liquid phase. This structure is used to define the advective term in
   * tracer equations.
   */

  cs_gwf_darcy_flux_t          *darcy;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var moisture_content
   * This quantity describes the level of saturation in a soil. This is a
   * constant value in case of a satured soil and variable one in case of a
   * unsaturated soil.
   */

  cs_property_t                *moisture_content;

  /*!
   * @}
   * @name Additional fields/arrays
   * @{
   *
   * \var pressure_head
   * Pointer to a \ref cs_field_t structure. Allocated only if the gravitation
   * effect is active. Location of this field depends on the discretization
   * scheme used to solve the Richards equation.
   * The pressure head is denoted by h, hydraulic head (unknowns solved in the
   * Richards eq.) is denoted by H and there are linked by:
   * h = H - gravity_potential
   */

  cs_field_t                   *pressure_head;

  /*!
   * @}
   */

} cs_gwf_saturated_single_phase_t;

/*! \struct cs_gwf_unsaturated_single_phase_t
 *
 * \brief Structure to handle the modelling of a single-phase flows in a porous
 *        media considered as saturated or not. Several simplifications can be
 *        be operated in this context. Only the liquid phase is taken into
 *        account.
 */

typedef struct {

    /*!
   * @name Equation
   * @{
   *
   * \var richards
   * The Richards equation is the governing equation which corresponds to the
   * mass conservation of water. "hydraulic_head" is the associated variable
   * "permeability" is the diffusion property related to this equation.
   */

  cs_equation_t                *richards;

  /*!
   * @}
   * @name Darcy (advection) field
   * @{
   *
   * \var darcy
   * Pointer to a \ref cs_gwf_darcy_flux_t structure. Darcy advective flux in
   * the liquid phase. This structure is used to define the advective term in
   * tracer equations.
   */

  cs_gwf_darcy_flux_t          *darcy;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var permeability
   * This quantity the product of the absolute permeability and the relative
   * permeability
   *
   * \var moisture_content
   * This quantity describes the level of saturation in a soil. This is a
   * constant value in case of a satured soil and variable one in case of a
   * unsaturated soil.
   *
   * \var soil_capacity
   * property attached to the unsteady term in the Richards equation
   */

  cs_property_t                *permeability;
  cs_property_t                *moisture_content;
  cs_property_t                *soil_capacity;

  /*!
   * @}
   * @name Additional fields/arrays
   * @{
   *
   * \var permeability_field
   * Pointer to a \ref cs_field_t structure. May be not allocated according to
   * the set of options. Store the value of the full permeability field in each
   * cell (k = k_abs * k_rel)
   *
   * \var moisture_field
   * Pointer to a \ref cs_field_t structure. Structure storing the value of the
   * moisture content in each cell. This is an optional structure i.e. it may
   * be set to NULL (for instance in case of a satured soil on the full
   * domain).
   *
   * \var capacity_field
   * Pointer to a \ref cs_field_t structure. Structure storing the value of the
   * soil capacity in each cell. This is an optional structure i.e. it is
   * set to NULL in case of a satured soil on the full domain.
   *
   * \var pressure_head
   * Pointer to a \ref cs_field_t structure. Allocated only if the gravitation
   * effect is active. Location of this field depends on the discretization
   * scheme used to solve the Richards equation.
   * The pressure head is denoted by h, hydraulic head (unknowns solved in the
   * Richards eq.) is denoted by H and there are linked by:
   * h = H - gravity_potential
   *
   * \var head_in_law
   * Array of values located at the same location as the hydraulic head solved
   * in the Richards equation. The values stored in this array are used to
   * update related quantities with law such as Van Genuchten/Mualen.
   */

  cs_field_t                   *permeability_field;
  cs_field_t                   *moisture_field;
  cs_field_t                   *capacity_field;
  cs_field_t                   *pressure_head;

  cs_real_t                    *head_in_law;

  /*!
   * @}
   */

} cs_gwf_unsaturated_single_phase_t;

/*! \struct cs_gwf_two_phase_t
 *
 * \brief Structure to handle the modelling of miscible or immiscible two-phase
 *        flows in a porous media.

 * The model follows what is depicted in "Finite volume approximation of a
 * diffusion-dissolution model and application to nuclear waste storage"
 * O. Angelini, C. Chavant, E. Chénier, R. Eymard and S. Granet in Mathematics
 * and Computers in Simulation (2011), 81 (10), pp. 2001--2017
 *
 * Main assumptions are:
 *   - No water in the gaseous phase
 *   - Incompressibility of the liquid phase
 *   - Hydrogen pressure is given by the "perfect gas" law in the gas phase and
 *     the Henry's law in the liquid phase
 *
 * The two primitive variables are the liquid and gas pressures with a specific
 * treatment in the saturated case to handle the gas pressure (cf. the cited
 * article or Angelini's PhD thesis)
 *
 * Notations are the following :
 * - Two phases: Liquid phase denoted by "l" and gaseous phase denoted by g
 * - Two components: water denoted by w and a gaseous component (let's say
 *   hydrogen) denoted by h. The gaseous component is present in the two phases
 *   whereas water is only considered in the liquid phase.
 *
 * The resulting linear algebraic system (one applies a linearization) is
 * defined as follows:
 *
 *                         liquid   gas
 * water conservation    | M_wl  | M_wg ||P_l|   | b_w |
 *                       |-------|------||---| = |-----|
 * hydrogen conservation | M_hl  | M_hg ||P_g|   | b_h |
 *
 * This is coupled system. Coupling terms are collected inside M_wg and M_hl
 */

typedef struct {

  /* Set of equations associated to this modelling */

  /*!
   * @name Equations and system of equations
   * @{
   *
   * \var wl_eq
   * Equation of conservation for the water component. Only the liquid phase is
   * considered. One assumes no water vapour in the gaseous phase. This
   * corresponds to the M_wl  block in the system of equations and to the b_w
   * right-hand side.
   */

  cs_equation_t                *wl_eq;

  /*! \var hg_eq
   * Equation of conservation for the (di)hydrogen. Hydrogen can be present in
   * the liquid or in the gaseous phase. This corresponds to the block (1,1) in
   * the system of equations, i.e. the M_hg block and the b_g right-hand side
   */

  cs_equation_t                *hg_eq;

  /*! \var wg_eqp
   * Parameters associated to the M_wg block i.e. the (0,1) block in the system
   * of equations. Water conservation w.r.t. the pressure in the gaseous phase.
   */

  cs_equation_param_t          *wg_eqp;

  /*! \var hl_eqp
   * Parameters associated to the (h,l) block i.e. the (1,0) block in the
   * system of equations. Conservation of the hydrogen w.r.t. the pressure in
   * the liquid phase.
   */

  cs_equation_param_t          *hl_eqp;

  /*! \var system
   * System of equations (wl_eq, hg_eq and the cross-term defined in the related
   * cs_equation_param_t structures) used for the coupled approach
   */

  cs_equation_system_t         *system;

  /*!
   * @}
   * @name Darcy (advection) fields
   * @{
   *
   * \var l_darcy
   * Pointer to a \ref cs_gwf_darcy_flux_t structure. Darcy advective flux in
   * the liquid phase
   *
   * \var g_darcy
   * Pointer to a \ref cs_gwf_darcy_flux_t structure. Darcy advective flux in
   * the gas phase
   */

  cs_gwf_darcy_flux_t          *l_darcy;
  cs_gwf_darcy_flux_t          *g_darcy;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var time_wl_pty
   * Property related to the unsteady term of the water conservation equation
   * w.r.t. the pressure in the liquid phase
   *
   * \var diff_wl_pty
   * Property related to the diffusion term of the water conservation equation
   * w.r.t. the pressure in the liquid phase
   *
   * \var time_wg_pty
   * Property related to the unsteady term of the water conservation equation
   * w.r.t. the pressure in the gaseous phase
   *
   * \var time_hg_pty
   * Property related to the unsteady term of the hydrogen conservation equation
   * w.r.t. the pressure in the gaseous phase.
   *
   * \var diff_hg_pty
   * Property related to the diffusion term of the hydrogen conservation
   * equation w.r.t. the pressure in the gaseous phase
   *
   * \var reac_hg_pty
   * Property related to the reaction term of the hydrogen conservation
   * equation w.r.t. the pressure in the gaseous phase
   *
   * \var time_hl_pty
   * Property related to the unsteady term of the hydrogen conservation equation
   * w.r.t. the pressure in the liquid phase.
   *
   * \var diff_hl_pty
   * Property related to the diffusion term of the hydrogen conservation
   * equation w.r.t. the pressure in the liquid phase
   *
   */

  cs_property_t                *time_wl_pty;
  cs_property_t                *diff_wl_pty;

  cs_property_t                *time_wg_pty;

  cs_property_t                *time_hg_pty;
  cs_property_t                *diff_hg_pty;
  cs_property_t                *reac_hg_pty;

  cs_property_t                *time_hl_pty;
  cs_property_t                *diff_hl_pty;

  /*!
   * @}
   * @name Additional fields
   * @{
   *
   * \var l_saturation
   *      Pointer to a \ref cs_field_t structure. Liquid saturation at cells.
   *      This quantity is denoted by \f$ S_l \f$ and is defined by the soil
   *      model
   *
   * \var c_pressure
   *      Pointer to a \ref cs_field_t structure named "capillarity_pressure".
   *      Capillarity pressure \f$ P_c = P_g - P_l \f$
   *
   * \var l_pressure
   *      Pointer to a \ref cs_field_t structure named "liquid_pressure".
   *      Pressure in the liquid phase is denoted by \f$ P_l \f$.
   *
   * \var g_pressure
   *      Pointer to a \ref cs_field_t structure named "gas_pressure".
   *      Pressure in the gas phase is denoted by \f$ P_g \f$.
   */

  cs_field_t                   *l_saturation;
  cs_field_t                   *c_pressure;
  cs_field_t                   *l_pressure;
  cs_field_t                   *g_pressure;

  /*!
   * @}
   * @name Additional arrays
   * @{
   *
   * \var time_wl_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the water conservation equation for the liquid
   *      phase. This array is linked to the \ref time_wl_pty (size = n_cells)
   *
   * \var diff_wl_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the water conservation equation. This array is linked to the
   *      \ref diff_wl_pty (size = n_cells)
   *
   * \var srct_wl_array
   *      Values in each cell of the source term corresponding to a quantity
   *      proportional to the time evolution of the liquid saturation
   *
   * \var time_wg_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the water conservation equation w.r.t. the pressure
   *      in the gaseous phase. This array is linked to the \ref time_wg_pty
   *      (size = n_cells)
   *
   * \var time_hg_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the hydrogen conservation equation w.r.t. the
   *      pressure in the gaseous phase. This array is linked to the \ref
   *      time_hg_pty (size = n_cells)
   *
   * \var diff_hg_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the hydrogen conservation equation w.r.t. to the pressure in
   *      the gaseous phase.  This array is linked to the \ref diff_hg_pty
   *      (size = n_cells)
   *
   * \var reac_hg_array
   *      Values in each cell of the coefficient appearing in the reaction
   *      term in the hydrogen conservation equation w.r.t. to the pressure in
   *      the gaseous phase.  This array is linked to the \ref reac_hg_pty
   *      (size = n_cells)
   *
   * \var srct_hg_array
   *      Values in each vertex of the source term corresponding to a quantity
   *      proportional to the time evolution of the unsteady term in the hg
   *      equation. This term is allocated only if the flag \ref
   *      CS_GWF_LIQUID_SATURATION_ON_SUBMESH is activated
   *
   * \var time_hl_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the hydrogen conservation equation w.r.t. the
   *      pressure in the liquid phase. This array is linked to the \ref
   *      time_hl_pty (size = n_cells)
   *
   * \var diff_hl_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the hydrogen conservation equation w.r.t. the pressure in the
   *      liquid phase.  This array is linked to the \ref diff_hl_pty (size
   *      = n_cells)
   *
   * \var l_rel_permeability
   *      Values in each cell of the relative permeability in the liquid phase.
   *      This quantity is used either in the water conservation or in the
   *      hydrogen conservation. This enables also to recover the (full)
   *      permeability in the liquid phase since
   *      permeability = abs_permeability * rel_l_permeability
   *      This quantity is defined by the soil model.
   *
   * \var g_rel_permeability
   *      Values in each cell of the relative permeability in the gas phase.
   *      This quantity is used either in the water conservation or in the
   *      hydrogen conservation. This enables also to recover the (full)
   *      permeability in the gas phase since
   *      permeability = abs_permeability * rel_l_permeability
   *      This quantity is defined by the soil model.
   *
   * \var l_capacity
   *      Values in each cell of the soil capacity defined as
   *      \f$ \frac{\partial S_l}{\partial P_c} \f$
   *      This quantity is defined by the soil model.
   *
   * \var capillarity_cell_pressure
   *      Values in each cell of the capillarity pressure. This quantity is the
   *      one used to update the variable related to a soil model such as the
   *      liquid and gaseous relative permeabilities or the liquid saturation.
   *
   * \var g_cell_pressure
   *      Values in each cell of the gas pressure. This quantity is the one
   *      used to update the unsteady/diffusion term coefficients in the
   *      conservation equation for the gas component. This quantity can also
   *      be useful for post-processing purposes.
   *
   * \var l_saturation_submesh
   *      Array storing the current values of the liquid saturation on a
   *      submesh. This submesh corresponds to the subdivision of the primal
   *      mesh by the dual mesh associated to each vertex (scanned this array
   *      with the c2v adjacency structure). This array is allocated only if
   *      the option CS_GWF_LIQUID_SATURATION_ON_SUBMESH is swicth on.
   *
   * \var l_saturation_submesh _pre
   *      Array storing the previous values of the liquid saturation on a
   *      submesh. This submesh corresponds to the subdivision of the primal
   *      mesh by the dual mesh associated to each vertex (scanned this array
   *      with the c2v adjacency structure). This array is allocated only if
   *      the option CS_GWF_LIQUID_SATURATION_ON_SUBMESH is swicth on.
   */

  cs_real_t                    *time_wl_array;
  cs_real_t                    *diff_wl_array;
  cs_real_t                    *srct_wl_array;

  cs_real_t                    *time_wg_array;

  cs_real_t                    *time_hg_array;
  cs_real_t                    *diff_hg_array;
  cs_real_t                    *reac_hg_array;
  cs_real_t                    *srct_hg_array;

  cs_real_t                    *time_hl_array;
  cs_real_t                    *diff_hl_array;

  cs_real_t                    *l_rel_permeability;
  cs_real_t                    *g_rel_permeability;

  /* These arrays are not always allocated. It depends on the numerical
     settings */

  cs_real_t                    *l_capacity;
  cs_real_t                    *capillarity_cell_pressure;
  cs_real_t                    *g_cell_pressure;
  cs_real_t                    *l_saturation_submesh;
  cs_real_t                    *l_saturation_submesh_pre;

  /*!
   * @}
   * @name Physical model parameters
   * @{
   *
   * \var l_mass_density
   *      Mass density in the liquid phase. With the model assumptions, this
   *      corresponds to the mass density of the main component in the liquid
   *      phase (e.g. water) in kg.m^-3
   *
   * \var l_viscosity
   *      Viscosity in the liquid phase (assumed to be constant) in Pa.s
   *
   * \var g_viscosity
   *      Viscosity in the gas phase (assumed to be constant) in Pa.s
   *
   * \var l_diffusivity_h
   *      Molecular diffusivity of the hydrogen in the liquid phase in m^2.s^-1
   *
   * \var w_molar_mass
   *      Molar mass of the main component in the liquid phase (e.g. water) in
   *      kg.mol^-1
   *
   * \var h_molar_mass
   *      Molar mass of the main component in the gas phase (e.g. hydrogen) in
   *      kg.mol^-1
   *
   * \var ref_temperature
   *      Reference temperature used in the "perfect gas" law (this is used
   *      when no thermal equation is solved). One expects a temperature in
   *      Kelvin.
   *
   * \var henry_constant
   *      Value of the Henry constant used in the Henry's law. Setting a very
   *      low value for this constant enables the model to degenerate into an
   *      immiscible model.
   */

  cs_real_t                     l_mass_density;
  cs_real_t                     l_viscosity;
  cs_real_t                     g_viscosity;
  cs_real_t                     l_diffusivity_h;
  cs_real_t                     w_molar_mass;
  cs_real_t                     h_molar_mass;
  cs_real_t                     ref_temperature;
  cs_real_t                     henry_constant;

  /*!
   * @}
   * @name Numerical parameters
   * @{
   *
   * \var use_coupled_solver
   * \brief When a model relies on several coupled equations, there are two
   *        main options to build and solve the system of equations. Either use
   *        a coupled solver (and thus build a coupled system) or use a
   *        segregated approach and an associated strategy to solve the
   *        sequence of equations and apply sub-iterations. The latter case
   *        (segregated solver) corresponds to the default choice.  true if a
   *        coupled solver is used. Otherwise a segregated solver is considered
   *
   * \var use_incremental_solver
   * \brief When a model includes non-linearities it can be useful to formulate
   *        the problem using increment and to iterate on the non-linear
   *        process (for instance whith a Picard or Anderson acceleration)
   *
   * \var use_properties_on_submesh
   * \brief Consider a submesh to define the liquid saturation
   *
   * \var nl_algo_type
   *      Type of algorithm to solve the non-linearities
   *
   * \var nl_algo_param
   *      Set of parameters to control the non-linear solver
   *
   * \var anderson_param
   *      Set of parameters to drive the Anderson acceleration (useful if the
   *      type of non-linear algorithm is set to the Anderson acceleration).
   *
   * \var nl_algo
   *      Structure used to manage the non-linearities
   */

  bool                           use_coupled_solver;
  bool                           use_incremental_solver;
  bool                           use_properties_on_submesh;

  cs_param_nl_algo_t             nl_algo_type;
  cs_iter_algo_param_t           nl_algo_param;
  cs_iter_algo_param_aa_t        anderson_param;

  cs_iter_algo_t                *nl_algo;

  /*!
   * @}
   */

} cs_gwf_two_phase_t;


/*! \struct cs_gwf_t
 *
 * \brief Main set of parameters/structures to manage the groundwater flow
 *        (GWF) module. This is an explicit definition of the structure
 *        \ref cs_gwf_t
 */

typedef struct {

  /*!
   * @name Metadata
   * @{
   *
   * \var model
   *      Model used to describe the behavior of the flow in the GWF module
   *      (system of equations related to the chosen physical modelling). See
   *      \ref cs_gwf_model_type_t for more details on each model
   *
   * \var flag
   *      Flag dedicated to general options to handle the GWF module
   *
   * \var post_flag
   *      Flag dedicated to the (automatic) post-processing of the GWF module
   */

  cs_gwf_model_type_t           model;
  cs_flag_t                     flag;
  cs_flag_t                     post_flag;

  /*!
   * @}
   * @name Properties
   * @{
   *
   * \var porosity
   *      Also called the saturated moisture content. This is a
   *      characterization of the portion of volume in a soil where the liquid
   *      (or also the gas) can be present. All models relies on this quantity.
   *
   * \var abs_permeability
   *      Absolute (or intrinsic) permeability which characterizes the behavior
   *      of a soil. According to the model of soil, this absolute permeability
   *      can be weigthed by a relative (scalar-valued) permeability. In the
   *      simplest case (saturated soil) the relative permeability is useless
   *      since this is equal to 1 (no weigth).
   */

  cs_property_t                *soil_porosity;
  cs_property_t                *abs_permeability;

  /*!
   * @}
   * @name Modelling context
   * @{
   *
   * \var model_context
   *      Pointer to a structure cast on-the-fly which depends on the choice of
   *      model (for instance single-phase or two-phase flows in porous media)
   */

  void                         *model_context;

  /*!
   * @}
   */

} cs_gwf_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in]      loc_flag   flag to define where the flux is defined
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_darcy_flux_t *
cs_gwf_darcy_flux_create(cs_flag_t         loc_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] p_darcy   pointer of pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_free(cs_gwf_darcy_flux_t  **p_darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in, out] darcy   pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_log(cs_gwf_darcy_flux_t    *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the definition of the advection field attached to a
 *         \ref cs_gwf_darcy_flux_t structure
 *
 * \param[in]       connect       pointer to a cs_cdo_connect_t structure
 * \param[in]       quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]       space_scheme  space discretization using this structure
 * \param[in, out]  darcy         pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_define(const cs_cdo_connect_t       *connect,
                         const cs_cdo_quantities_t    *quant,
                         cs_param_space_scheme_t       space_scheme,
                         cs_gwf_darcy_flux_t          *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the advection field/arrays related to the Darcy flux.
 *
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      eq        pointer to the equation related to this Darcy flux
 * \param[in]      cur2prev  true or false
 * \param[in, out] darcy     pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_update(const cs_real_t              t_eval,
                         const cs_equation_t         *eq,
                         bool                         cur2prev,
                         cs_gwf_darcy_flux_t         *darcy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate the balance by zone (relying on the splitting arising from
 *         the boundary settings) for the advection field attached to a \ref
 *         cs_gwf_darcy_flux_t structure
 *
 * \param[in]       connect       pointer to a cs_cdo_connect_t structure
 * \param[in]       quant         pointer to a cs_cdo_quantities_t structure
 * \param[in]       eqp           pointer to the set of equation parameters
 * \param[in, out]  darcy         pointer to the darcy structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_darcy_flux_balance(const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant,
                          const cs_equation_param_t    *eqp,
                          cs_gwf_darcy_flux_t          *darcy);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PRIV_H__ */
