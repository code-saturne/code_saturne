#ifndef __CS_GWF_HYDRAULIC_MODEL_H__
#define __CS_GWF_HYDRAULIC_MODEL_H__

/*============================================================================
 * Structures related to the hydraulic models available in the GWF module
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_equation_system.h"
#include "cs_gwf_priv.h"
#include "cs_param_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* ++++++++++++++++++++++++++++++++++++++++++++++ */
/* Saturated single-phase flows in a porous media */
/* ++++++++++++++++++++++++++++++++++++++++++++++ */

/*! \struct cs_gwf_sspf_t
 *
 * \brief Structure to handle the modelling of a single-phase flows in a porous
 *        media considered as saturated.
 *
 * Several simplifications are operated in this context. Only the liquid phase
 * is taken into account. Several properties are constants.
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
   * constant value in case of a satured soil.
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

} cs_gwf_sspf_t;

/* ++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Unsaturated single-phase flows in a porous media */
/* ++++++++++++++++++++++++++++++++++++++++++++++++ */

/*! \struct cs_gwf_uspf_t
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
   * moisture content in each cell.
   *
   * \var capacity_field
   * Pointer to a \ref cs_field_t structure. Structure storing the value of the
   * soil capacity in each cell.
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

} cs_gwf_uspf_t;

/* +++++++++++++++++++++++++++++++++ */
/* Two-phase flows in a porous media */
/* +++++++++++++++++++++++++++++++++ */

/*! \struct cs_gwf_tpf_t
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
 * The two primitive variables are the capillarity and liquid pressures with a
 * specific treatment in the saturated case (cf. the cited article or
 * Angelini's PhD thesis)
 *
 * Notations are the following :
 * - Two phases: Liquid phase denoted by "l" and gaseous phase denoted by "g"
 * - indice "c" refers to the capillarity pressure
 * - Two components: water denoted by "w" and a gaseous component (let's say
 *   hydrogen) denoted by "h". The gaseous component is present in the two
 *   phases whereas water is only considered in the liquid phase.
 *
 * The resulting linear algebraic system (one applies a linearization) is
 * defined as follows:
 *
 *                              cap.    liq
 * water mass conservation    | M_00  | M_01 ||P_c|   | b_w |
 *                            |-------|------||---| = |-----|
 * hydrogen mass conservation | M_10  | M_11 ||P_l|   | b_h |
 *
 * This is a coupled system. Coupling terms are collected inside M_01 and M_10
 */

typedef struct {

  /* Set of equations associated to this modelling */

  /*!
   * @name Equations and system of equations
   * @{
   *
   * \var w_eq
   * Equation of mass conservation for water. Only the liquid phase is
   * considered. One assumes no water vapour in the gaseous phase. This
   * corresponds to the M_00 and M_01 blocks in the system of equations and to
   * the b_w right-hand side.
   */

  cs_equation_t                *w_eq;

  /*! \var h_eq
   * Equation of mass conservation for (di)hydrogen. Hydrogen can be present in
   * the liquid or in the gaseous phase. This corresponds to the M_10 and M_11
   * blocks in the system of equations along with the b_h right-hand side
   */

  cs_equation_t                *h_eq;

  /*! \var b01_w_eqp
   * Parameters associated to the (0,1) block in the system of equations. Water
   * conservation w.r.t. the capillarity pressure.
   */

  cs_equation_param_t          *b01_w_eqp;

  /*! \var b10_h_eqp
   * Parameters associated to the (1,0) block in the system of equations.
   * Conservation of the hydrogen w.r.t. the capillarity pressure.
   */

  cs_equation_param_t          *b10_h_eqp;

  /*! \var system
   * System of equations (w_eq, h_eq and the cross-terms defined in the related
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
   *
   * \var t_darcy
   * Pointer to a \ref cs_gwf_darcy_flux_t structure. Darcy advective flux for
   * the total flux (linear combination of the liquid/gas Darcy flux)
   */

  cs_gwf_darcy_flux_t          *l_darcy;
  cs_gwf_darcy_flux_t          *g_darcy;
  cs_gwf_darcy_flux_t          *t_darcy;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var time_wc_pty
   * Property related to the unsteady term of the water conservation equation
   * w.r.t. the capillarity pressure
   *
   * \var diff_wl_pty
   * Property related to the diffusion term of the water conservation equation
   * w.r.t. the pressure in the liquid phase
   *
   * \var time_hc_pty
   * Property related to the unsteady term of the hydrogen conservation equation
   * w.r.t. the capillarity pressure
   *
   * \var diff_hc_pty
   * Property related to the diffusion term of the hydrogen conservation
   * equation w.r.t. the capillarity pressure
   *
   * \var time_hl_pty
   * Property related to the unsteady term of the hydrogen conservation equation
   * w.r.t. the pressure in the liquid phase.
   *
   * \var diff_hl_pty
   * Property related to the diffusion term of the hydrogen conservation
   * equation w.r.t. the pressure in the liquid phase
   *
   * \var reac_h_pty
   * Property related to the reaction term of the hydrogen conservation
   * equation w.r.t. the pressure in the gas phase. Only used when a segregated
   * solver is considered.
   *
   * \var diff_g_pty
   * Property used in the definition of the Darcy flux in the gas phase
   */

  cs_property_t                *time_wc_pty;
  cs_property_t                *diff_wl_pty;

  cs_property_t                *time_hc_pty;
  cs_property_t                *diff_hc_pty;

  cs_property_t                *time_hl_pty;
  cs_property_t                *diff_hl_pty;

  cs_property_t                *reac_h_pty;

  cs_property_t                *diff_g_pty;

  /*!
   * @}
   * @name Additional fields
   * @{
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
   *
   * \var l_saturation
   *      Pointer to a \ref cs_field_t structure. Liquid saturation at cells.
   *      This quantity is denoted by \f$ S_l \f$ and is defined by the soil
   *      model
   */

  cs_field_t                   *c_pressure;
  cs_field_t                   *l_pressure;
  cs_field_t                   *g_pressure;
  cs_field_t                   *l_saturation;

  /*!
   * @}
   * @name Additional arrays
   * @{
   *
   * \var time_wc_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the water conservation equation associated to the
   *      capillarity pressure as variable. This array is linked to the \ref
   *      time_wc_pty (size = n_cells or c2v->idx[n_cells] if the definition
   *      relies on a submesh)
   *
   * \var diff_wl_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the water conservation equation. This array is linked to the
   *      \ref diff_wl_pty (size = n_cells)
   *
   * \var time_hc_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the hydrogen conservation equation associated to the
   *      capillarity pressure. This array is linked to the \ref time_hc_pty
   *      (size = n_cells or c2v->idx[n_cells] if the definition relies on a
   *      submesh)
   *
   * \var diff_hc_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the hydrogen conservation equation associated to the
   *      capillarity pressure This array is linked to the \ref diff_hc_pty
   *      (size = n_cells). This term may be useless according to the numerical
   *      options.
   *
   * \var time_hl_array
   *      Values in each cell of the coefficient appearing in front of the
   *      unsteady term in the hydrogen conservation equation w.r.t. the
   *      pressure in the liquid phase. This array is linked to the \ref
   *      time_hl_pty (size = n_cells or c2v->idx[n_cells] if the definition
   *      relies on a submesh)
   *
   * \var diff_hl_array
   *      Values in each cell of the coefficient appearing in the diffusion
   *      term in the hydrogen conservation equation w.r.t. the pressure in the
   *      liquid phase.  This array is linked to the \ref diff_hl_pty (size
   *      = n_cells)
   *
   * \var reac_h_array
   *      Values of the reaction coefficient appearing in the diffusion term in
   *      the hydrogen conservation equation w.r.t. the pressure in the liquid
   *      phase.  This array is linked to the \ref diff_hl_pty (size = n_cells
   *      or c2v->idx[n_cells] if the definition relies on a submesh)
   *
   * \var srct_w_array
   *      Values of the source terms for the water conservation equation. Only
   *      used if a segregated solver is considered. Size = n_cells or
   *      c2v->idx[n_cells] if the definition relies on a submesh
   *
   * \var srct_h_array
   *      Values of the source terms for the hydrogen conservation equation.
   *      Only used if a segregated solver is considered. Size = n_cells or
   *      c2v->idx[n_cells] if the definition relies on a submesh
   *
   * \var reac_h_array
   *      Values of the reaction coefficient appearing in the hydrogen
   *      conservation equation. Only used if a segregated solver is
   *      considered. Size = n_cells or c2v->idx[n_cells] if the definition
   *      relies on a submesh
   *
   * \var diff_g_array
   *      Values of the diffusion coefficient used to compute the Darcy flux
   *      in the gaz phase. Size = n_cells
   */

  cs_real_t                    *time_wc_array;
  cs_real_t                    *diff_wl_array;

  cs_real_t                    *time_hc_array;
  cs_real_t                    *diff_hc_array;

  cs_real_t                    *time_hl_array;
  cs_real_t                    *diff_hl_array;

  cs_real_t                    *srct_w_array;
  cs_real_t                    *srct_h_array;
  cs_real_t                    *reac_h_array;

  cs_real_t                    *diff_g_array;

  /*
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
   * \var c_pressure_cells
   *      Values in each cell of the capillarity pressure. This quantity is the
   *      one used to update the variable related to a soil model such as the
   *      liquid and gaseous relative permeabilities or the liquid saturation.
   *
   * \var g_pressure_cells
   *      Values in each cell of the gas pressure. This quantity is the one
   *      used to update the unsteady/diffusion term coefficients in the
   *      conservation equation for the gas component. This quantity can also
   *      be useful for post-processing purposes.
   *
   * \var l_capacity
   *      Values in each cell of the soil capacity defined as
   *      \f$ \frac{\partial S_l}{\partial P_c} \f$
   *      This quantity is defined by the soil model.
   *
   * \var l_saturation_submesh
   *      Array storing the current values of the liquid saturation on a
   *      submesh. This submesh corresponds to the subdivision of the primal
   *      mesh by the dual mesh associated to each vertex (scanned this array
   *      with the c2v adjacency structure). This array is allocated only if
   *      the option CS_GWF_LIQUID_SATURATION_ON_SUBMESH is switch on.
   */

  cs_real_t                    *l_rel_permeability;
  cs_real_t                    *g_rel_permeability;

  /* These arrays are not always allocated. It depends on the numerical
     settings */

  cs_real_t                    *c_pressure_cells;
  cs_real_t                    *g_pressure_cells;
  cs_real_t                    *l_capacity;
  cs_real_t                    *l_saturation_submesh;

  /*!
   * @}
   * @name Physical model parameters
   * @{
   *
   * \var is_miscible
   *      true or false. Some simplifications can be done if immiscible.
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

  bool                          is_miscible;

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
   * \var use_definition_on_submesh
   * \brief Consider a submesh to define the liquid saturation
   *
   * \var use_diffusion_view_for_darcy
   * \brief Use a diffusion term for the discretization of the Darcy terms in
   *        the conservation equation for the mass of hydrogen. The default
   *        option is to consuder an advection term since it should be more
   *        robust as upwinding technique can be used.
   *
   * \var nl_algo_type
   *      Type of algorithm to solve the non-linearities
   *
   * \var nl_algo_cvg
   *      Set of parameters to drive the convergence of the non-linear solver
   *
   * \var nl_relax_factor
   *      Value of the relaxation factor in the non-linear algorithm. A classical
   *      choice is between 0.70 and 0.95
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
  bool                           use_definition_on_submesh;
  bool                           use_diffusion_view_for_darcy;

  cs_param_nl_algo_t             nl_algo_type;
  cs_param_sles_cvg_t            nl_algo_cvg;
  cs_real_t                      nl_relax_factor;
  cs_iter_algo_param_aa_t        anderson_param;

  cs_iter_algo_t                *nl_algo;

  /*!
   * @}
   */

} cs_gwf_tpf_t;


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
   *
   * \var verbosity
   *      level of information printed
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

  int                           verbosity;
  cs_gwf_model_type_t           model;
  cs_flag_t                     flag;
  cs_flag_t                     post_flag;

  /*!
   * @}
   * @name Properties
   * @{
   *
   * \var soil_porosity
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
   * @name Hydraulic modelling context
   * @{
   *
   * \var model_context
   *      Pointer to a structure cast on-the-fly which depends on the choice of
   *      the hydraulic model (for instance single-phase or two-phase flows in
   *      porous media)
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

END_C_DECLS

#endif /* __CS_GWF_HYDRAULIC_MODEL_H__ */
