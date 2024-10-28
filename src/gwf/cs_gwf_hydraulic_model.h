#ifndef __CS_GWF_HYDRAULIC_MODEL_H__
#define __CS_GWF_HYDRAULIC_MODEL_H__

/*============================================================================
 * Structures related to the hydraulic models available in the GWF module
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
   * update related quantities with law such as Van Genuchten/Mualem.
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

/*!
 * \enum cs_gwf_tpf_solver_type_t
 * \brief Type of solver considered for a two-phase flow model
 *
 * \var CS_GWF_TPF_SOLVER_PLPC_COUPLED
 *      (Pc, Pl) is the couple of main unknowns. Pl for the conservation of the
 *      mass of water and Pc for the conservation of the mass of component
 *      mainly present in the gas phase (H2 for instance). A fully coupled
 *      approach solver is considered to solve the system of equations.
 *
 * \var CS_GWF_TPF_SOLVER_PLPC_COUPLED_INCR
 *      (dPc, dPl) is the couple of main unknowns. dPc and dPl represents the
 *      increment on the capillarity pressure and liquid pressure
 *      respectively. The two equations solved are the conservation of the mass
 *      of water and the conservation of the mass of component mainly present
 *      in the gas phase (H2 for instance). A fully coupled approach solver is
 *      considered to solve the system of equations along with a modified
 *      Picard algorithm in the spirit of Celia.
 */

typedef enum {

  CS_GWF_TPF_SOLVER_PLPC_COUPLED,
  CS_GWF_TPF_SOLVER_PLPC_COUPLED_INCR,

  CS_GWF_TPF_N_SOLVERS

} cs_gwf_tpf_solver_type_t;

/*!
 * \enum cs_gwf_tpf_approx_type_t
 * \brief Type of approximation used for the computation of unsteady or
 *        diffusion properties
 */

typedef enum {

  CS_GWF_TPF_APPROX_PC_CELL_AVERAGE,
  CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE,
  CS_GWF_TPF_APPROX_PC_EDGE_AVERAGE,
  CS_GWF_TPF_APPROX_PC_VERTEX_AVERAGE,

  CS_GWF_TPF_APPROX_VERTEX_SUBCELL,

  CS_GWF_TPF_N_APPROX

} cs_gwf_tpf_approx_type_t;

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
 *     the Henry's law in the liquid phase (when a miscible model is used)
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
   */

  cs_gwf_darcy_flux_t          *l_darcy;
  cs_gwf_darcy_flux_t          *g_darcy;

  /*!
   * @}
   * @name Properties related to the model
   * @{
   *
   * \var krl_pty
   * Property related to the relative permeability in the liquid phase
   *
   * \var krg_pty
   * Property related to the relative permeability in the gas phase
   *
   * \var lsat_pty
   * Property related to the liquid saturation
   *
   * \var lcap_pty
   * Property related to the liquid capacity (derivative of the liquid
   * saturation w.r.t. the capillarity pressure)
   */

  cs_property_t                *krl_pty;
  cs_property_t                *krg_pty;
  cs_property_t                *lsat_pty;
  cs_property_t                *lcap_pty;

  /* Properties associated to a discrete term in the system of equations
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
   * \var diff_g_pty
   * Property used in the definition of the Darcy flux in the gas phase
   */

  cs_property_t                *time_wc_pty;
  cs_property_t                *diff_wl_pty;

  cs_property_t                *time_hc_pty;
  cs_property_t                *diff_hc_pty;

  cs_property_t                *time_hl_pty;
  cs_property_t                *diff_hl_pty;

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
   *
   * \var l_rho_h
   *      Pointer to a \ref cs_field_t structure. Property related to the mass
   *      density in the liquid phase of the component mainly present in the
   *      gas phase (e.g. H2). The component present in the liquid comes from
   *      the dissolution process (ruled by the Henry law)
   *
   * \var g_rho_h
   *      Pointer to a \ref cs_field_t structure. Property related to the mass
   *      density in the gas phase of the component mainly present in the gas
   *      phase (e.g. H2). The component present in the gas phase is ruled by
   *      the perfect gas law.
   */

  cs_field_t                   *c_pressure;
  cs_field_t                   *l_pressure;
  cs_field_t                   *g_pressure;
  cs_field_t                   *l_saturation;
  cs_field_t                   *l_rho_h;
  cs_field_t                   *g_rho_h;

  /*!
   * @}
   * @name Additional arrays
   * @{
   *
   * \var srct_w_array
   *      Values of the source terms for the water conservation equation.
   *      Size = n_vertices
   *
   * \var srct_h_array
   *      Values of the source terms for the hydrogen conservation equation.
   *      Size = n_vertrices
   *
   * \var lsat_pre_array
   *      Values of the liquid saturation at the previous time step when an
   *      approximation \ref CS_GWF_TPF_APPROX_VERTEX_SUBCELL has been chosen
   *      Size = c2v index size
   */

  cs_real_t                    *srct_w_array;
  cs_real_t                    *srct_h_array;
  cs_real_t                    *lsat_pre_array;

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
  cs_real_t                     h_molar_mass;
  cs_real_t                     ref_temperature;
  cs_real_t                     henry_constant;

  /*!
   * @}
   * @name Numerical parameters
   * @{
   *
   * \var solver_type
   * \brief Type of solver considered to solve the system of equations (choice
   *        of main unknowns and strategy of resolution (coupled/segregated))
   *
   * \var enforce_pg_positivity
   * \brief Apply a rescaling on the gas pressure to avoid a negative value
   *
   * \var approx_type
   * \brief type of approximation used for the computation of diffusion,
   *        unsteady coefficients
   *
   * \var cell_weight
   * \brief weight associated to the cell value w.r.t. to the values at
   *        vertices when a \ref CS_GWF_TPF_APPROX_PC_CELL_VERTEX_AVERAGE
   *        choice is set for 'approx_type'. If the value is < 0 or > 1, then
   *        the default settings is kept.
   *
   * \var upwind_weight
   * \brief weight associated to the upwind value w.r.t. to the centered value
   *        There is no effect when \ref CS_GWF_TPF_APPROX_VERTEX_SUBCELL is
   *        chosen to approximate coefficients. If the value is < 0 or > 1,
   *        then the default settings is kept.
   *
   * \var nl_algo_type
   *      Type of algorithm to solve the non-linearities
   *
   * \var nl_algo_verbosity
   *      Level of verbosity for the algorithm solving the non-linearities
   *
   * \var nl_cvg_param
   *      Set of parameters to drive the convergence of the non-linear solver
   *
   * \var anderson_param
   *      Set of parameters to drive the Anderson acceleration (useful if the
   *      type of non-linear algorithm is set to the Anderson acceleration).
   *
   * \var nl_algo
   *      Structure used to manage the non-linearities
   */

  cs_gwf_tpf_solver_type_t       solver_type;
  bool                           enforce_pg_positivity;
  cs_gwf_tpf_approx_type_t       approx_type;
  double                         cell_weight;
  double                         upwind_weight;

  cs_param_nl_algo_t             nl_algo_type;
  int                            nl_algo_verbosity;
  cs_param_convergence_t         nl_cvg_param;
  cs_iter_algo_param_aac_t       anderson_param;

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
