#ifndef __CS_NAVSTO_PARAM_H__
#define __CS_NAVSTO_PARAM_H__

/*============================================================================
 * Routines to handle cs_navsto_param_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_boundary.h"
#include "cs_equation_param.h"
#include "cs_math.h"
#include "cs_physical_constants.h"
#include "cs_sles.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Manage the naming of properties, variables and equations related to the
 * Navier-Stokes module */

#define CS_NAVSTO_LAMINAR_VISCOSITY  "laminar_viscosity"
#define CS_NAVSTO_STREAM_EQNAME      "streamfunction_eq"

/*!
 * @name Flags specifying numerical options
 * @{
 */

/* Value = 1 */
#define CS_NAVSTO_FLAG_STEADY            (1 <<  0) /*!< Steady-state */

/*!
 * @}
 * @name Flag specifying predefined post-processing
 *
 * \brief w denotes the vorticity * vector and u the velocity vector, k is the
 *        kinetic energy defined by k := 1/2 * u \cdot u
 *
 * @{
 */

/* Value =   1 */
#define CS_NAVSTO_POST_VELOCITY_DIVERGENCE (1 <<  0) /*!< div(u) */

/* Value =   2 */
#define CS_NAVSTO_POST_KINETIC_ENERGY      (1 <<  1) /*!< k := rho/2 u\cdot u */

/* Value =   4 */
#define CS_NAVSTO_POST_VORTICITY           (1 <<  2) /*!< w = curl(u) */

/* Value =   8 */
#define CS_NAVSTO_POST_VELOCITY_GRADIENT   (1 <<  3)

/* Value =  16 */
#define CS_NAVSTO_POST_STREAM_FUNCTION     (1 <<  4) /*!< -Lap(Psi) = w_z */

/* Value =  32 */
#define CS_NAVSTO_POST_HELICITY            (1 <<  5) /*!< u \cdot w */

/* Value =  64 */
#define CS_NAVSTO_POST_ENSTROPHY           (1 <<  6) /*!< w \cdot w */

/*!
 * @}
 */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef cs_flag_t  cs_navsto_param_model_t;

/*! \enum cs_navsto_param_model_bit_t
 *  \brief Bit values for physical modelling related to the Navier-Stokes system
 *  of equations
 *
 * \var CS_NAVSTO_MODEL_STOKES
 * Stokes equations (mass and momentum) with the classical choice of variables
 * i.e. velocity and pressure. Mass density is assumed to be constant.
 *
 * \var CS_NAVSTO_MODEL_OSEEN
 * Like the incompressible Navier-Stokes equations (mass and momentum) but with
 * a velocity field which is given. Thus the advection term in the momentum
 * equation is linear. Unknowns: velocity and pressure. Mass density is assumed
 * to be constant.
 *
 * \var CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
 * Navier-Stokes equations: mass and momentum with a constant mass density
 *
 * \var CS_NAVSTO_MODEL_GRAVITY_EFFECTS
 * Take into account the gravity effects (add a constant source term equal to
 * rho*vect(g))
 *
 * \var CS_NAVSTO_MODEL_CORIOLIS_EFFECTS
 * Take into account the Coriolis effects (add a source term)
 *
 * \var CS_NAVSTO_MODEL_BOUSSINESQ
 * Gravity effects are taken into account as well as the effect of small
 * variation of temperatures.
 * The gradient of temperature is assumed to have a small norm and the mass
 * density variates in a small range. In this case, an additional equation
 * related to the temperature is considered and momentum source term is added.
 */

typedef enum {

  /* Main modelling for the dynamic
     ------------------------------ */

  CS_NAVSTO_MODEL_STOKES                          = 1<<0, /* =   1 */
  CS_NAVSTO_MODEL_OSEEN                           = 1<<1, /* =   2 */
  CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES    = 1<<2, /* =   4 */

  /* Additional modelling bits
     ------------------------- */

  CS_NAVSTO_MODEL_GRAVITY_EFFECTS                 = 1<<3, /* =   8 */
  CS_NAVSTO_MODEL_CORIOLIS_EFFECTS                = 1<<4, /* =  16 */
  CS_NAVSTO_MODEL_BOUSSINESQ                      = 1<<5  /* =  32 */

} cs_navsto_param_model_bit_t;

/*! \enum cs_navsto_sles_t
 *
 *  \brief High-level information about the way of settings the SLES for solving
 *  the Navier-Stokes system. When the system is treated as a saddle-point
 *  problem (monolithic approach in what follows), then one uses these
 *  notations: A_{00} is the upper-left block and A_{11} (should be 0 but the
 *  preconditioner may have entries for the approximation of the inverse of the
 *  Schur complement).
 *
 * \var CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK
 * Associated keyword: "no_block"
 *
 * Use the same mechanism as for a stand-alone equation. In this case, the
 * setting relies on the function \ref cs_equation_set_sles and the different
 * options for solving a linear system such as the choice of the iterative
 * solver or the choice of the preconditioner or the type of residual
 * normalization
 *
 *
 * \var CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG
 * Associated keyword: "block_amg_cg"
 *
 * The Navier-Stokes system of equations is solved using a multigrid on each
 * diagonal block as a preconditioner and applying a conjugate gradient as
 * solver. Use this strategy when the saddle-point problem has been reformulated
 * into a "classical" linear system. For instance when a Uzawa or an Artificial
 * Compressibility coupling algorithm is used. (i.e. with the parameter
 * \ref CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY or
 * \ref CS_NAVSTO_COUPLING_UZAWA is set as coupling algorithm). This option is
 * only available with the support to the PETSc library up to now.
 *
 *
 * \var CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK
 * Associated keyword: "additive_gmres"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm) The Navier-Stokes
 * system of equations is solved an additive preconditioner (block diagonal
 * matrix where the block 00 is A_{00}) and the block 11 is set to the identity.
 * Preconditioner/solver for the block 00 is set using the momentum equation.
 * This option is only available with the support to the PETSc library up to now.
 *
 *
 * \var CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK
 * Associated keyword: "multiplicative_gmres"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm) The Navier-Stokes
 * system of equations is solved a multiplicative preconditioner (block diagonal
 * matrix where the block 00 is A_{00}) and the block 11 is set to the identity.
 * Block 01 is also considered in the block preconditioner.
 * Preconditioner/solver for the block 00 is set using the momentum equation.
 * This option is only available with the support to the PETSc library up to now.
 *
 *
 * \var CS_NAVSTO_SLES_DIAG_SCHUR_GMRES
 * Associated keyword: "diag_schur_gmres"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm). The
 * Navier-Stokes system of equations is solved using a block diagonal
 * preconditioner where the block 00 is A_{00} preconditioned with one multigrid
 * iteration and the block 11 is an approximation of the Schur complement
 * preconditioned with one multigrid iteration. The main iterative solver is a
 * flexible GMRES. This option is only available with the support to the PETSc
 * library up to now.
 *
 *
 * \var CS_NAVSTO_SLES_UPPER_SCHUR_GMRES
 * Associated keyword: "upper_schur_gmres"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm). The
 * Navier-Stokes system of equations is solved using a upper triangular block
 * preconditioner where the block 00 is A_{00} preconditioned with one multigrid
 * iteration and the block 11 is an approximation of the Schur complement
 * preconditioned with a minres. The main iterative solver is a flexible
 * GMRES. This option is only available with the support to the PETSc
 * library up to now.
 *
 * \var CS_NAVSTO_SLES_GKB
 * Associated keyword: "gkb"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm). The
 * Navier-Stokes system of equations is solved using a Golub-Kahan
 * bi-diagonalization. One assumes that the saddle-point system is symmetric.
 * By default, the block A_{00} may be augmented (this is not the default
 * choice) and is solved with a conjugate gradient algorithm preconditioned
 * with a multigrid. The residual is computed in the energy norm. This option is
 * only available with the support to the PETSc library up to now.
 *
 * * \var CS_NAVSTO_SLES_GKB_GMRES
 * Associated keyword: "gkb_gmres"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm). The
 * Navier-Stokes system of equations is solved using a Golub-Kahan
 * bi-diagonalization (GKB) as preconditioner of a flexible GMRES solver. The
 * GKB algorithm is solved with a reduced tolerance as well as the CG+Multigrid
 * used as an inner solver in the GKB algorithm. One assumes that the
 * saddle-point system is symmetric. The residual for the GKB part is computed
 * in the energy norm. This option is only available with the support to the
 * PETSc library up to now.
 *
 * \var CS_NAVSTO_SLES_GKB_SATURNE
 * Associated keyword: "gkb_saturne"
 *
 * Available choice when a monolithic approach is used (i.e. with the parameter
 * CS_NAVSTO_COUPLING_MONOLITHIC is set as coupling algorithm). The
 * Navier-Stokes system of equations is solved using a Golub-Kahan
 * bi-diagonalization.
 * By default, the block A_{00} may be augmented (this is not the default
 * choice) and is solved with the SLES settings given to the momentum equation
 * A conjugate gradient algorithm preconditioned
 * with a multigrid for Stokes for instance.
 * The residual is computed in the energy norm.
 *
 * \var CS_NAVSTO_SLES_MUMPS
 * Associated keyword: "mumps"
 *
 * Direct solver to solve systems arising from the discretization of the
 * Navier-Stokes equations
 *
 * \var CS_NAVSTO_SLES_UZAWA_AL
 * Associated keyword: "uzawa_al"
 * Resolution using an uzawa algorithm with an Augmented Lagrangian approach
 */

typedef enum {

  CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK,
  CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG,
  CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK,
  CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK,
  CS_NAVSTO_SLES_DIAG_SCHUR_GMRES,
  CS_NAVSTO_SLES_UPPER_SCHUR_GMRES,
  CS_NAVSTO_SLES_GKB,
  CS_NAVSTO_SLES_GKB_GMRES,
  CS_NAVSTO_SLES_GKB_SATURNE,
  CS_NAVSTO_SLES_MUMPS,
  CS_NAVSTO_SLES_UZAWA_AL,

  CS_NAVSTO_SLES_N_TYPES

} cs_navsto_sles_t;


/*! \struct cs_navsto_param_sles_t
 *  \brief Structure storing the parameters for solving the Navier-Stokes system
 */

typedef struct {

  /*! \var strategy
   * Choice of strategy for solving the Navier--Stokes system
   */
  cs_navsto_sles_t              strategy;

  /*! \var algo_tolerance
   *  Tolerance at which the Oseen/Stokes system is resolved (apply to the
   *  residual of the coupling algorithm chosen to solve the Navier--Stokes
   *  system)
   */
  cs_real_t                     algo_tolerance;

  /*! \var algo_n_max_iter
   * Maximal number of iterations of the coupling algorithm.
   */
  int                           algo_n_max_iter;

  /*! \var picard_tolerance
   *  Tolerance at which the Picard algorithm is resolved. One handles the
   *  non-linearity arising from the advection term with the algorithm.
   */
  cs_real_t                     picard_tolerance;

  /*! \var picard_n_max_iter
   * Maximal number of iterations for the Picard algorithm used to handle
   * the non-linearity arising from the advection term.
   */
  int                           picard_n_max_iter;

} cs_navsto_param_sles_t;


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

  /*! \var verbosity
   * Level of display of the information related to the Navier-Stokes system
   */
  int                         verbosity;

  /*! \var post_flag
   * Flag storing which predefined post-processing has to be done
   */
  cs_flag_t                   post_flag;

  /*!
   * @name Physical modelling
   * Which equations to solve ?  Properties and their related fields are
   * allocated according to the choice of model for Navier-Stokes
   * @{
   */

  /*! \var model
   * Modelling related to the Navier-Stokes system of equations
   */
  cs_navsto_param_model_t     model;

  /*! \var reference_pressure
   *  Value of the reference pressure p0 (used for rescaling or during update
   *  of physical quantities). By default: 0.
   */

  cs_real_t                   reference_pressure;

  /*! \var phys_constants
   * Main physical constants (gravity vector and coriolis source term). This
   * structure is shared with the legacy part.
   */
  cs_physical_constants_t    *phys_constants;

  /*! \var density
   *  Density of the fluid, pointer to \ref cs_property_t used in several
   *  terms in the Navier-Stokes equations
   */

  cs_property_t              *density;

  /*! \var lami_viscosity
   *  Laminar viscosity, pointer to \ref cs_property_t associated to the
   *  diffusion term for the momentum equation
   */

  cs_property_t              *lami_viscosity;

  /*!
   * @}
   * @name Turbulence modelling
   * Set of parameters to handle turbulence modelling.
   * @{
   */

  /*! \var turbulence
   * Main set of parameters to handle turbulence modelling. This
   * structure is shared with the legacy part.
   */

  cs_turb_model_t            *turbulence;

  /*! \var rans_modelling
   * Main set of parameters to handle RANS modelling. This
   * structure is shared with the legacy part.
   * RANS means Reynolds Average Navier-Stokes
   */

  cs_turb_rans_model_t       *rans_modelling;


  /*!
   * @name Numerical options
   * Set of numerical options to build the linear system and how to solve it
   * @{
   */

  /*! \var option_flag
   * Flag storing high-level option related to the Navier-Stokes system
   */
  cs_flag_t                     option_flag;

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

  /*! \var dof_reduction_mode
   *  How are defined the Degrees of freedom
   */
  cs_param_dof_reduction_t      dof_reduction_mode;

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

  /*! \var adv_form
   *  Type of formulation for the advection term
   *
   *  \var adv_scheme
   *  Type of scheme for the advection term
   */
  cs_param_advection_form_t     adv_form;
  cs_param_advection_scheme_t   adv_scheme;

  /*! \var qtype
   *  A \ref cs_quadrature_type_t indicating the type of quadrature to use in
   *  all routines involving quadratures
   */
  cs_quadrature_type_t          qtype;


  /*! \var sles_param
   * Set of choices to control the resolution of the Navier--Stokes system
   */
  cs_navsto_param_sles_t        sles_param;

  /*!
   * @}
   * @name Initial conditions (IC)
   *
   * Set of parameters used to take into account the initial condition on the
   * pressure and/or the velocity.
   * CAUTION: so far, there is no check if the different IC are compatible
   * with the boundary conditions for instance
   * @{
   */

  /*! \var velocity_ic_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside a \ref cs_equation_param_t
   *  structure dedicated to the momentum equation.
   *
   * \var n_velocity_ic_defs
   *  Number of initial conditions associated to the velocity
   *
   * \var velocity_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  velocity.
   *  The code does not check if the resulting initial velocity satisfies the
   *  divergence constraint.
   */

  bool         velocity_ic_is_owner;
  int          n_velocity_ic_defs;
  cs_xdef_t  **velocity_ic_defs;

  /*! \var pressure_ic_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside a dedicated \ref cs_equation_param_t
   *
   * \var n_pressure_ic_defs
   *  Number of initial conditions associated to the pressure
   *
   * \var pressure_ic_defs
   *  Pointers to the definitions of the initial conditions associated to the
   *  pressure.
   *  In order to force a zero-mean pressure, the code can compute the average
   *  of the resulting pressure and subtract it
   */

  bool         pressure_ic_is_owner;
  int          n_pressure_ic_defs;
  cs_xdef_t  **pressure_ic_defs;

  /*! @}
   * @name Boundary conditions (BC)
   *
   * Set of parameters used to take into account the boundary conditions on the
   * pressure and/or the velocity.
   * @{
   */

  /* \var boundaries
   * Pointer to a \ref cs_boundary_t structure shared with the domain
   */
  const cs_boundary_t   *boundaries;

  /*! \var velocity_bc_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside a dedicated \ref cs_equation_param_t
   *  Most of the time this should be false since an equation is associated to
   *  to the resolution of the velocity field (the momentum equation).
   *
   * \var n_velocity_bc_defs
   * Number of definitions related to the settings of the boundary conditions
   * for the velocity field.
   *
   * \var velocity_bc_defs
   * Array of pointers to the definition of boundary conditions for the velocity
   * field
   */

  bool         velocity_bc_is_owner;
  int          n_velocity_bc_defs;
  cs_xdef_t  **velocity_bc_defs;

  /*! \var pressure_bc_is_owner
   *  True if the definitions are stored inside this structure, otherwise
   *  the definitions are stored inside a dedicated \ref cs_equation_param_t
   *  if an equation solves the pressure field.
   *
   * \var n_pressure_bc_defs
   *  Number of boundary conditions associated to the pressure field.
   *
   * \var pressure_bc_defs
   *  Pointers to the definitions of the boundary conditions associated to the
   *  pressure field.
   */

  bool         pressure_bc_is_owner;
  int          n_pressure_bc_defs;
  cs_xdef_t  **pressure_bc_defs;

  /*! @} */

} cs_navsto_param_t;

/*! \struct cs_navsto_algo_info_t
 *  \brief Set of information related to the convergence of the iterative
 *         algorithm (Picard or Uzawa for instance)
 *
 * \var cvg
 * convergence or divergence status
 *
 * \var res
 * value of the residual for the iterative algorithm
 *
 * \var n_algo_iter
 * number of iterations for the algorithm (outer iterations)
 *
 * \var n_inner_iter
 * cumulated number of inner iterations (sum over the outer iterations)
 *
 * \var last_inner_iter
 * last number of iterations for the inner solver
 */

typedef struct {

  cs_sles_convergence_state_t      cvg;
  double                           res;
  int                              n_algo_iter;
  int                              n_inner_iter;
  int                              last_inner_iter;

} cs_navsto_algo_info_t;

/*! \enum cs_navsto_key_t
 *  \brief List of available keys for setting the parameters of the
 *         Navier-Stokes system
 *
 * \var CS_NSKEY_ADVECTION_FORMULATION
 * Set the type of formulation for the advection term, for example in the  Oseen
 * problem . cf. cs_param.h
 *
 * \var CS_NSKEY_ADVECTION_SCHEME
 * Set the type of scheme for the advection term, for example in the  Oseen
 * problem . cf. cs_param.h
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
 * Set the maximal number of iteration for solving the coupled system.
 *
 * \var CS_NSKEY_MAX_PICARD_ITER
 * Set the maximal number of Picard iterations for solving the non-linearity
 * arising from the advection form
 *
 * \var CS_NSKEY_QUADRATURE
 * Set the type to use in all routines involving quadrature (similar to \ref
 * CS_EQKEY_BC_QUADRATURE)
 *
 * \var CS_NSKEY_PICARD_TOLERANCE
 * Tolerance at which the non-linearity arising from the advection term is
 * resolved
 *
 * \var CS_NSKEY_RESIDUAL_TOLERANCE
 * Tolerance at which the Oseen or Stokes system is resolved (apply to the
 * residual of the coupling algorithm chosen to solve the Navier--Stokes system)
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

  CS_NSKEY_ADVECTION_FORMULATION,
  CS_NSKEY_ADVECTION_SCHEME,
  CS_NSKEY_DOF_REDUCTION,
  CS_NSKEY_GD_SCALE_COEF,
  CS_NSKEY_MAX_ALGO_ITER,
  CS_NSKEY_MAX_PICARD_ITER,
  CS_NSKEY_QUADRATURE,
  CS_NSKEY_PICARD_TOLERANCE,
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
 * \brief  Initialize a cs_navsto_algo_info_t structure
 *
 * \param[in]   ns_info   pointer to a cs_navsto_algo_info_t to initialize
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_navsto_algo_info_init(cs_navsto_algo_info_t   *ns_info)
{
  if (ns_info == NULL)
    return;

  ns_info->cvg = CS_SLES_ITERATING;
  ns_info->res = cs_math_big_r;
  ns_info->n_algo_iter = 0;
  ns_info->n_inner_iter = 0;
  ns_info->last_inner_iter = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print header before dumping information gathered in the structure
 *         cs_navsto_algo_info_t
 *
 * \param[in]  algo_name     name of the algorithm
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_navsto_algo_info_header(const char   *algo_name)
{
  assert(algo_name != NULL);
  cs_log_printf(CS_LOG_DEFAULT,
                "%8s.It  -- Algo.Res   Inner    Cumul  ||div(u)||\n", algo_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print header before dumping information gathered in the structure
 *         cs_navsto_algo_info_t
 *
 * \param[in]  algo_name     name of the algorithm
 * \param[in]  ns_info       cs_navsto_algo_info_t structure
 * \param[in]  div_l2        l2 norm of the divergence
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_navsto_algo_info_printf(const char                    *algo_name,
                           const cs_navsto_algo_info_t    ns_info,
                           double                         div_l2)
{
  assert(algo_name != NULL);
  if (ns_info.n_algo_iter == 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  "%8s.It%02d-- %9s  %5d  %6d  %6.4e\n",
                  algo_name, ns_info.n_algo_iter, " ",
                  ns_info.last_inner_iter, ns_info.n_inner_iter, div_l2);
  else
    cs_log_printf(CS_LOG_DEFAULT,
                  "%8s.It%02d-- %5.3e  %5d  %6d  %6.4e\n",
                  algo_name, ns_info.n_algo_iter, ns_info.res,
                  ns_info.last_inner_iter, ns_info.n_inner_iter, div_l2);
}

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
cs_navsto_param_is_steady(const cs_navsto_param_t       *nsp)
{
  if (nsp == NULL)
    return true;

  if (nsp->option_flag & CS_NAVSTO_FLAG_STEADY)
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
 * \param[in]  algo_coupling  algorithm used for solving the NS system
 * \param[in]  option_flag    additional high-level numerical options
 * \param[in]  post_flag      predefined post-processings
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t             *boundaries,
                       cs_navsto_param_model_t          model,
                       cs_navsto_param_coupling_t       algo_coupling,
                       cs_flag_t                        option_flag,
                       cs_flag_t                        post_flag);

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
 * \brief  Set the value to consider for the reference pressure
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]  pref      value of the reference pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_reference_pressure(cs_navsto_param_t    *nsp,
                                 cs_real_t             pref);

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
 * \brief  Define the initial condition for the velocity unknowns.
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
 * \brief  Define the initial condition for the pressure unknowns.
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
 * \brief  Add the definition of boundary conditions related to a fixed wall
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_fixed_walls(cs_navsto_param_t    *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to a symmetry
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_symmetries(cs_navsto_param_t    *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to outlets
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_outlets(cs_navsto_param_t    *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the pressure field on a boundary using a uniform value.
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      value     value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_pressure_bc_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for a sliding wall boundary using a
 *         uniform value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_wall_by_value(cs_navsto_param_t    *nsp,
                                     const char           *z_name,
                                     cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a uniform
 *         value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_inlet_by_value(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an analytical
 *         function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_inlet_by_analytic(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_analytic_func_t   *ana,
                                         void                 *input);

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
                                   bool                  is_owner,
                                   cs_lnum_t            *index);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a advection field for the Oseen problem
 *
 * \param[in, out]    nsp        pointer to a \ref cs_navsto_param_t
 * \param[in, out]    adv_fld    pointer to a \ref cs_adv_field_t
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_oseen_field(cs_navsto_param_t   *nsp,
                          cs_adv_field_t      *adv_fld);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_PARAM_H__ */
