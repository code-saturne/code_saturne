#ifndef __CS_NAVSTO_PARAM_H__
#define __CS_NAVSTO_PARAM_H__

/*============================================================================
 * Routines to handle cs_navsto_param_t structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_navsto_param_model_t
 *  \brief Modelling related to the Navier-Stokes system of equations
 *
 * \var CS_NAVSTO_MODEL_STOKES
 * Stokes equations (mass and momentum) with the classical choice of variables
 * i.e. velocity and pressure
 *
 * \var CS_NAVSTO_MODEL_OSEEN
 * Like the incompressible Navier-Stokes equations (mass and momentum) but with
 * a velocity field which is given. Thus the advection term in the momentum
 * equation is linear. Unknowns: velocity and pressure
 *
 * \var CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
 * Navier-Stokes equations: mass and momentum with a constant mass density
 *
 * \var CS_NAVSTO_MODEL_BOUSSINESQ_NAVIER_STOKES
 * Navier-Stokes equations: mass and momentum with a constant mass density
 * The gradient of temperature is assumed to have a small norm and the mass
 * density variates in a small range. In this case, an additional equation
 * related to the energy is considered.
 */

typedef enum {

  CS_NAVSTO_MODEL_STOKES,
  CS_NAVSTO_MODEL_OSEEN,
  CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
  CS_NAVSTO_MODEL_BOUSSINESQ_NAVIER_STOKES,

  CS_NAVSTO_N_MODELS

} cs_navsto_param_model_t;

/*! \enum cs_navsto_param_sles_t
 *  \brief High-level information about the Way of settings the SLES
 *  for solving the Navier-Stokes system
 *
 * \var CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK
 * Use the same mechanism as for stand-alone equation relying on the function
 * \ref cs_equation_set_sles
 *
 * \var CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG
 * The Navier-Stokes system of equations is solved using a multigrid on each
 * diagonal block as a preconditioner and applying a conjugate gradient as
 * solver. Use this strategy when the saddle-point problem has been reformulated
 * into a "classical" linear system. For instance when a Uzawa or an Artificial
 * Compressibility coupling algorithm is used. This option is only available
 * with the support to the PETSc library up to now.
 *
 * \var CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK
 * The Navier-Stokes system of equations is solved an additive preconditioner
 * (block diagonal matrix where the block 00 is A_{00} preconditionned by one
 * multigrid iteration and the block 11 is set to the identity. This option
 * is only available with the support to the PETSc library up to now.
 * Available choice when a monolithic approach is used.
 *
 * \var CS_NAVSTO_SLES_DIAG_SCHUR_GMRES
 *
 * The Navier-Stokes system of equations is solved using a block diagonal
 * preconditioner where the block 00 is A_{00} preconditioned with one multigrid
 * iteration and the block 11 is an approximation of the Schur complement
 * preconditionned with one multigrid iteration. The main iterative solver is a
 * flexible GMRES. Available choice when a monolithic approach is used.
 *
 * \var CS_NAVSTO_SLES_UPPER_SCHUR_GMRES
 *
 * The Navier-Stokes system of equations is solved using a upper triangular
 * block preconditioner where the block 00 is A_{00} preconditioned with one
 * multigrid iteration and the block 11 is an approximation of the Schur
 * complement preconditionned with a minres. The main iterative solver is a
 * flexible GMRES. Available choice when a monolithic approach is used.
 */

typedef enum {

  CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK,
  CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG,
  CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK,
  CS_NAVSTO_SLES_DIAG_SCHUR_GMRES,
  CS_NAVSTO_SLES_UPPER_SCHUR_GMRES,

  CS_NAVSTO_SLES_N_TYPES

} cs_navsto_param_sles_t;

/*! \enum cs_navsto_param_time_state_t
 *  \brief Status of the time for the Navier-Stokes system of equations
 *
 * \var CS_NAVSTO_TIME_STATE_UNSTEADY
 * The Navier-Stokes system of equations is time-dependent
 *
 * \var CS_NAVSTO_TIME_STATE_FULL_STEADY
 * The Navier-Stokes system of equations is solved without taking into account
 * the time effect
 *
 * \var CS_NAVSTO_TIME_STATE_LIMIT_STEADY
 * The Navier-Stokes system of equations is solved as a limit of a unsteady
 * process
 */

typedef enum {

  CS_NAVSTO_TIME_STATE_FULL_STEADY,
  CS_NAVSTO_TIME_STATE_LIMIT_STEADY,
  CS_NAVSTO_TIME_STATE_UNSTEADY,

  CS_NAVSTO_N_TIME_STATES

} cs_navsto_param_time_state_t;

/*! \enum cs_navsto_param_coupling_t
 *  \brief Choice of algorithm for solving the system
 *
 * \var CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY
 * The system is solved using an artificial compressibility algorithm.
 * One vectorial equation is solved followed by a pressure update.
 *
 * \var CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP
 * The system is solved using an artificial compressibility algorithm with a
 * Vector Penalty Projection splitting.
 * Two vectorial equations are solved: a momentum-like one and another one
 * involving a grad-div operator.
 *
 * \var CS_NAVSTO_COUPLING_MONOLITHIC
 * The system is treated as a "monolithic" matrix
 *
 * \var CS_NAVSTO_COUPLING_PROJECTION
 * The system is solved using an incremental projection algorithm
 *
 * \var CS_NAVSTO_COUPLING_UZAWA
 * The system is solved without decoupling the equations using a Uzawa algorithm
 * and an Augmented Lagrangian approach inside each sub-iteration.
 */

typedef enum {

  CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY,
  CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP,
  CS_NAVSTO_COUPLING_MONOLITHIC,
  CS_NAVSTO_COUPLING_PROJECTION,
  CS_NAVSTO_COUPLING_UZAWA,

  CS_NAVSTO_N_COUPLINGS

} cs_navsto_param_coupling_t;

/*! \struct cs_navsto_param_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         Navier-Stokes system
 */

typedef struct {

  /* \var boundaries
   * Pointer to a \ref cs_boundary_t structure shared with the domain
   */
  const cs_boundary_t          *boundaries;

  /*! \var verbosity
   * Level of display of the information related to the Navier-Stokes system
   */
  int                           verbosity;

  /*! \var dof_reduction_mode
   *  How are defined the Degrees of freedom
   */
  cs_param_dof_reduction_t      dof_reduction_mode;

  /*! \var time_scheme
   * Discretization scheme for time
   *
   * \var theta
   * Value of the parameter for the time scheme when a theta-scheme is used
   *
   */
  cs_param_time_scheme_t        time_scheme;
  cs_real_t                     theta;

  /*! \var space_scheme
   * Discretization scheme for space
   */
  cs_param_space_scheme_t       space_scheme;

  /*! \var model
   * Modelling related to the Navier-Stokes system of equations
   */
  cs_navsto_param_model_t       model;

  /*!
   * \var has_gravity
   * Take into account the gravity effect: true or false
   *
   * \var gravity
   * Vector related to the gravity effect
   */
  bool                          has_gravity;
  cs_real_3_t                   gravity;

  /*! \var time_state
   * Status of the time for the Navier-Stokes system of equations
   */
  cs_navsto_param_time_state_t  time_state;

  /*! \var sles_strategy
   * Choice of strategy for solving the SLES system
   */
  cs_navsto_param_sles_t        sles_strategy;

  /*! \var coupling
   * Choice of algorithm for solving the system
   */
  cs_navsto_param_coupling_t    coupling;

  /*! \var gd_scale_coef
   *  Default value to set the scaling of the grad-div term when an
   *  artificial compressibility algorithm or an Uzawa-Augmented Lagrangian
   *  method is used
   */
  cs_real_t                     gd_scale_coef;

  /*! \var qtype
   *  A \ref cs_quadrature_type_t indicating the type of quadrature to use in
   *  all routines involving quadratures
   */
  cs_quadrature_type_t          qtype;

  /*! \var residual_tolerance
   *  Tolerance at which the Navier--Stokes is resolved (apply to the residual
   * of the coupling algorithm chosen to solve the Navier--Stokes system)
   */
  cs_real_t                     residual_tolerance;

  /*! \var max_algo_iter
   * Maximal number of iteration of the coupling algorithm. Not useful for a
   * monolithic approach. In this case, only the maximal number of iterations
   * for the iterative solver is taken into account
   */
  int                           max_algo_iter;

  /*!
   * @}
   * @name Physical properties
   * Set of properties: properties and their related fields are allocated
   * according to the choice of model for Navier-Stokes
   * @{
   */

  /*! \var density
   *  Density of the fluid, pointer to \ref cs_property_t used in several
   *  terms in the Navier-Stokes equations
   */

  cs_property_t      *density;

  /*! \var lami_viscosity
   *  Laminar viscosity, pointer to \ref cs_property_t associated to the
   *  diffusion term for the momentum equation
   */

  cs_property_t      *lami_viscosity;

  /*!
   * @}
   * @name Initial conditions(IC)
   *
   * Set of parameters used to take into account the initial condition on the
   * pressure and/or the velocity.
   * CAUTION: so far, there is no check if the different IC are compatible
   * with the boundary conditions for instance
   * @{
   */

  /*! \var velocity_ic_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside the a \ref cs_equation_param_t
   *  structure dedicated to the momentum equation.
   */
  bool  velocity_ic_is_owner;

  /*! \var n_velocity_ic_defs
   *  Number of initial conditions associated to the velocity
   */
  int   n_velocity_ic_defs;

  /*! \var velocity_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  velocity.
   *  The code does not check if the resulting initial velocity satisfies the
   *  divergence constraint.
   */
  cs_xdef_t  **velocity_ic_defs;

  /*! \var pressure_ic_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside a dedicated \ref cs_equation_param_t
   */
  bool  pressure_ic_is_owner;

  /*! \var n_pressure_ic_defs
   *  Number of initial conditions associated to the pressure
   */
  int   n_pressure_ic_defs;

  /*! \var pressure_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  pressure.
   *  In order to force a zero-mean pressure, the code can compute the average
   *  of the resulting pressure and subtract it
   */
  cs_xdef_t  **pressure_ic_defs;

  /*! @} */

} cs_navsto_param_t;

/*! \enum cs_navsto_key_t
 *  \brief List of available keys for setting the parameters of the
 *         Navier-Stokes system
 *
 * \var CS_NSKEY_DOF_REDUCTION
 * Set how the DoFs are defined (similar to \ref CS_EQKEY_DOF_REDUCTION)
 * Enable to set this type of DoFs definition for all related equations
 *
 * \var CS_NSKEY_GD_SCALE_COEF
 * Set the scaling of the grad-div term when an artificial compressibility
 * algorithm or an Uzawa - Augmented Lagrangian method is used
 *
 * \var CS_NSKEY_MAX_ALGO_ITER
 * Set the maximal number of iteration of the coupling algorithm. Not useful
 * for a monolithic approach. In this case, only the maximal number of
 * iterations for the iterative solver is taken into account
 *
 * \var CS_NSKEY_QUADRATURE
 * Set the type to use in all routines involving quadrature (similar to \ref
 * CS_EQKEY_BC_QUADRATURE)
 *
 * \var CS_NSKEY_RESIDUAL_TOLERANCE
 * Tolerance at which the Navier--Stokes is resolved (apply to the residual
 * of the coupling algorithm chosen to solve the Navier--Stokes system)
 *
 * \var CS_NSKEY_SLES_STRATEGY
 * Strategy for solving the SLES arising from the discretization of the
 * Navier-Stokes system
 *
 * \var CS_NSKEY_SPACE_SCHEME
 * Numerical scheme for the space discretization. Available choices are:
 * - "cdo_fb"  for CDO face-based scheme
 *
 * \var CS_NSKEY_TIME_SCHEME
 * Numerical scheme for the time discretization
 *
 * \var CS_NSKEY_TIME_THETA
 * Set the value of theta. Only useful if CS_NSKEY_TIME_SCHEME is set to
 * "theta_scheme"
 * - Example: "0.75" (keyval must be between 0 and 1)
 *
 * \var CS_NSKEY_VERBOSITY
 * Set the level of details for the specific part related to the Navier-Stokes
 * system
 */

typedef enum {

  CS_NSKEY_DOF_REDUCTION,
  CS_NSKEY_GD_SCALE_COEF,
  CS_NSKEY_MAX_ALGO_ITER,
  CS_NSKEY_QUADRATURE,
  CS_NSKEY_RESIDUAL_TOLERANCE,
  CS_NSKEY_SLES_STRATEGY,
  CS_NSKEY_SPACE_SCHEME,
  CS_NSKEY_TIME_SCHEME,
  CS_NSKEY_TIME_THETA,
  CS_NSKEY_VERBOSITY,

  CS_NSKEY_N_KEYS

} cs_navsto_key_t;

/*============================================================================
 * Inline static public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask \ref cs_navsto_param_t structure if the settings correspond to
 *         a steady computation
 *
 * \param[in]  nsp     pointer to a \ref cs_navsto_param_t structure
*
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_navsto_param_is_steady(cs_navsto_param_t       *nsp)
{
  if (nsp == NULL)
    return true;

  if (nsp->time_state == CS_NAVSTO_TIME_STATE_FULL_STEADY)
    return true;
  else
    return false;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in]  boundaries     pointer to a cs_boundary_t structure
 * \param[in]  model          model related to the NS system to solve
 * \param[in]  time_state     state of the time for the NS equations
 * \param[in]  algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t              *boundaries,
                       cs_navsto_param_model_t           model,
                       cs_navsto_param_time_state_t      time_state,
                       cs_navsto_param_coupling_t        algo_coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_navsto_param_t
 *         structure
 *
 * \param[in, out] nsp      pointer to a \ref cs_navsto_param_t structure to set
 * \param[in]      key      key related to the member of eq to set
 * \param[in]      keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set(cs_navsto_param_t    *nsp,
                    cs_navsto_key_t       key,
                    const char           *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the numerical settings defined for the Navier-Stokes system
 *         to an equation related to this system.
 *
 * \param[in]       nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in, out]  eqp    pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_transfer(const cs_navsto_param_t    *nsp,
                         cs_equation_param_t        *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of the coupling algorithm
 *
 * \param[in]     coupling    a \ref cs_navsto_param_coupling_t
 *
 * \return the name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_coupling_name(cs_navsto_param_coupling_t  coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_analytic(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_analytic_func_t   *ana,
                                      void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      val       pointer to the value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_val(cs_navsto_param_t    *nsp,
                                 const char           *z_name,
                                 cs_real_t            *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an array
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                           (true or false)
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_array(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_flag_t             loc,
                                   cs_real_t            *array,
                                   _Bool                 is_owner,
                                   cs_lnum_t            *index);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_PARAM_H__ */
