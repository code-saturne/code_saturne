#ifndef __CS_NAVSTO_PARAM_H__
#define __CS_NAVSTO_PARAM_H__

/*============================================================================
 * Functions to handle cs_navsto_param_t structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_boundary.h"
#include "cdo/cs_cdo_turbulence.h"
#include "cdo/cs_equation_param.h"
#include "cdo/cs_iter_algo.h"
#include "base/cs_math.h"
#include "alge/cs_param_sles.h"
#include "base/cs_physical_constants.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Manage the naming of properties, variables and equations related to the
 * Navier-Stokes module
 */

#define CS_NAVSTO_STREAM_EQNAME      "streamfunction_eq"

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef cs_flag_t  cs_navsto_param_model_flag_t;
typedef cs_flag_t  cs_navsto_param_post_flag_t;

/*! \enum cs_navsto_param_model_t
 *  \brief Describe the system of equations related to the Navier-Stokes to be
 *  solved
 *
 * \var CS_NAVSTO_MODEL_STOKES
 * Stokes equations (mass and momentum) with the classical choice of variables
 * i.e. velocity and pressure. Mass density is assumed to be constant.
 *
 * \var CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES
 * Navier-Stokes equations: mass and momentum with a constant mass density
 * Mass equation is equivalent to an incompressibility constraint
 *
 * \var CS_NAVSTO_MODEL_OSEEN
 * Like the incompressible Navier-Stokes equations (mass and momentum) but with
 * a velocity field which is given. Thus the advection term in the momentum
 * equation is linear. Unknowns: velocity and pressure. Mass density is assumed
 * to be constant. The advection field is set using
 * \ref cs_navsto_add_oseen_field
 *
 */

typedef enum {

  CS_NAVSTO_MODEL_STOKES,
  CS_NAVSTO_MODEL_OSEEN,
  CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,

  CS_NAVSTO_N_MODELS

} cs_navsto_param_model_t;

/*! \enum cs_navsto_param_model_bit_t
 *  \brief Bit values for additional physical modelling related to the
 *  Navier-Stokes system of equations
 *
 * \var CS_NAVSTO_MODEL_STEADY
 * There is no time dependency in the model
 *
 * \var CS_NAVSTO_MODEL_GRAVITY_EFFECTS
 * Take into account the gravity effects (add a constant source term equal to
 * rho*vect(g))
 *
 * \var CS_NAVSTO_MODEL_CORIOLIS_EFFECTS
 * Take into account the Coriolis effects (add a source term)
 *
 * \var CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER
 * An additional equation is created involving the thermal equation.
 * The advection field is automatically set as the mass flux. By default,
 * the temperature is solved in Celsius and homogeneous Neumann boundary
 * conditions (no flux) are set. This option is not compatible with the option
 * \ref CS_NAVSTO_MODEL_BOUSSINESQ
 *
 * \var CS_NAVSTO_MODEL_BOUSSINESQ
 * An additional equation is created involving the thermal equation. The
 * advection field is automatically set as the mass flux. By default, the
 * temperature is solved in Celsius and homogeneous Neumann boundary conditions
 * (no flux) are set.  The variation of mass density around a reference value is
 * related to the variation of temperature w.r.t. a reference temperature and it
 * plays a role in the gravity effects (rho.vect(g)) The gradient of temperature
 * is assumed to have a small norm and the mass density variates in a small
 * range. In this case, an additional momentum source term is added.
 *
 * \var CS_NAVSTO_MODEL_WITH_SOLIDIFICATION
 * The Navier-Stokes is modified to take into account solidification process.
 * A boussinesq term is added as well as a head loss term derived from a Darcy
 * approximation. Please see the \ref cs_solidification_model_t structure and
 * related structures for more details
 *
 */

typedef enum {

  CS_NAVSTO_MODEL_STEADY                          = 1<<0, /* =   1 */
  CS_NAVSTO_MODEL_GRAVITY_EFFECTS                 = 1<<1, /* =   2 */
  CS_NAVSTO_MODEL_CORIOLIS_EFFECTS                = 1<<2, /* =   4 */
  CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER          = 1<<3, /* =   8 */
  CS_NAVSTO_MODEL_BOUSSINESQ                      = 1<<4, /* =  16 */
  CS_NAVSTO_MODEL_WITH_SOLIDIFICATION             = 1<<5  /* =  32 */

} cs_navsto_param_model_bit_t;

/*! \enum cs_navsto_param_post_bit_t
 *  \brief Bit values for additional generic postprocessing related to the
 *  Navier-Stokes module. In what follows, w denotes the vorticity vector, u
 *  the velocity vector and k the kinetic energy defined by
 *  \f$ \frac{\rho u}{2} \cdot u \f$
 *
 * \var CS_NAVSTO_POST_VELOCITY_DIVERGENCE
 * Compute div(u) and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_KINETIC_ENERGY
 * Compute \f$ rho/2 u\cdot u \f$ and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_VORTICITY
 * Compute w = curl(u) and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_VELOCITY_GRADIENT
 * Compute the tensor grad(u) and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_STREAM_FUNCTION
 * Add an equation to compute the stream function and associate a field to this
 * quantity. This equation is a scalar valued diffusion equation:
 * -Lap(Psi) = w_z
 * This is only useful for a 2D computation (assuming that the z-azis is the
 * extruded one, i.e. the flow is in the x-y plane
 *
 * \var CS_NAVSTO_POST_HELICITY
 * Compute \f$ h = u \cdot w \f$ and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_ENSTROPHY
 * Compute \f$ w \cdot w \f$ and associate a field to this quantity
 *
 * \var CS_NAVSTO_POST_MASS_DENSITY
 * Compute the mass density in each cell and monitor the evolution of the
 * integral of the mass in the full computational domain. If a Boussinesq
 * approximation is set, then the mass density used in the buoyancy term is
 * considered (not the reference value assumed to be constant).
 *
 * \var CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE
 * Compute the balance of the mass flux for each cell. When the mass density is
 * constant, this quantity is only a scaling (by the mass density) of the
 * quantity computed with the flag CS_NAVSTO_POST_VELOCITY_DIVERGENCE.
 *
 * \var CS_NAVSTO_POST_PRESSURE_GRADIENT
 * Compute the presure gradient in each cell.
 */

typedef enum {

  CS_NAVSTO_POST_VELOCITY_DIVERGENCE      = 1<< 0, /* =   1 */
  CS_NAVSTO_POST_KINETIC_ENERGY           = 1<< 1, /* =   2 */
  CS_NAVSTO_POST_VORTICITY                = 1<< 2, /* =   4 */
  CS_NAVSTO_POST_VELOCITY_GRADIENT        = 1<< 3, /* =   8 */
  CS_NAVSTO_POST_STREAM_FUNCTION          = 1<< 4, /* =  16 */
  CS_NAVSTO_POST_HELICITY                 = 1<< 5, /* =  32 */
  CS_NAVSTO_POST_ENSTROPHY                = 1<< 6, /* =  64 */
  CS_NAVSTO_POST_MASS_DENSITY             = 1<< 7, /* = 128 */
  CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE   = 1<< 8, /* = 256 */
  CS_NAVSTO_POST_PRESSURE_GRADIENT        = 1<< 9, /* = 512 */

} cs_navsto_param_post_bit_t;

typedef struct {

  /*! \var verbosity
   *  Level of information printed
   */

  int                           verbosity;


} cs_navsto_param_sles_t;

/*! \enum cs_navsto_param_coupling_t
 *  \brief Choice of algorithm for solving the system
 *
 * \var CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY
 * The system is solved using an artificial compressibility algorithm.
 * One vectorial equation is solved followed by a pressure update.
 *
 * \var CS_NAVSTO_COUPLING_MONOLITHIC
 * The system is treated as a "monolithic" matrix
 *
 * \var CS_NAVSTO_COUPLING_PROJECTION
 * The system is solved using an incremental projection algorithm.
 * The potential is based on cell-based scheme
 */

typedef enum {

  CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY,
  CS_NAVSTO_COUPLING_MONOLITHIC,
  CS_NAVSTO_COUPLING_PROJECTION,

  CS_NAVSTO_N_COUPLINGS

} cs_navsto_param_coupling_t;

/*! \struct cs_navsto_param_boussinesq_t
 *  \brief Structure storing the parameters related to the Boussinesq source
 *         term in the momentum equation
 */

typedef struct {

  cs_real_t    beta;      /* Dilatation coefficient */
  cs_real_t    var0;      /* Reference value of the variable */

  /* Array of values of the variable (for instance the temperature). This is a
   * shared pointer. The lifecycle of this array is not managed by this
   * structure.
   */

  const cs_real_t   *var;

} cs_navsto_param_boussinesq_t;

/*! \struct cs_navsto_param_t
 *  \brief Structure storing the parameters related to the resolution of the
 *         Navier-Stokes system
 */

typedef struct {

  /*!
   * @name Physical modelling
   * Which equations to solve ?  Properties and their related fields are
   * allocated according to the choice of model for Navier-Stokes
   * @{
   */

  /*! \var model
   * Modelling related to the Navier-Stokes system of equations
   */

  cs_navsto_param_model_t        model;

  /*! \var model_flag
   * Flag storing high-level option related to the Navier-Stokes system
   */

  cs_navsto_param_model_flag_t   model_flag;

  /*! \var turbulence
   *  Structure storing all information needed to set the turbulence modelling
   */

  cs_turbulence_param_t         *turbulence;

  /*!
   * @}
   * @name Properties and fields related to the Navier-Stokes module
   * @{
   */

  /*! \var phys_constants
   * Main physical constants (gravity vector and Coriolis source term). This
   * structure is shared with the legacy part.
   */

  cs_physical_constants_t       *phys_constants;

  /*! \var mass_density
   * Mass_density of the fluid, pointer to \ref cs_property_t used in several
   * terms in the Navier-Stokes equations
   */

  cs_property_t                 *mass_density;

  /*! \var tot_viscosity
   *  Laminar viscosity + if needed the turbulent viscosity
   *  Pointer to \ref cs_property_t associated to the
   *  diffusion term for the momentum equation
   */

  cs_property_t                *tot_viscosity;

  /*! \var lam_viscosity
   *  Laminar viscosity
   */

  cs_property_t                *lam_viscosity;

  /*!
   * @}
   * @name Numerical options
   * Set of numerical options to build the linear system and how to solve it
   * @{
   */

  /*! \var dof_reduction_mode
   *  How are defined the Degrees of freedom
   */

  cs_param_dof_reduction_t       dof_reduction_mode;

  /*! \var coupling
   * Choice of algorithm for solving the system
   */

  cs_navsto_param_coupling_t     coupling;

  /*! \var space_scheme
   * Discretization scheme for space
   */

  cs_param_space_scheme_t        space_scheme;

  /* Boussinesq approximation:
   *
   * Take into account buoyancy terms (variation of mass density w.r.t. the
   * variation of a field (for instance the temperature but can be also a
   * concentration as in segregation model in the solidification module)
   *
   * \var n_boussinesq_terms
   * Number of contributions to the buoyancy source term in the Boussinesq
   * approximation
   *
   * \var boussinesq_param
   * Structure storing elements used to compute the Boussinesq approximation
   */

  int                                 n_boussinesq_terms;
  cs_navsto_param_boussinesq_t       *boussinesq_param;

  /*!
   * @}
   * @name Non-linear algorithm
   * Set of parameters to drive the resolution of the non-linearity arising from
   * the Navier--Stokes system
   * @{
   */

  /*! \var nl_algo_type
   *  Type of algorithm used to tackle the non-linearity arising from the
   *  advection term
   */

  cs_param_nl_algo_t            nl_algo_type;

  /*! \var nl_cvg_param
   *  Structure storing several tolerances and metadata to drive the
   *  convergence of the non-linear iterative algorithm used to solve te
   *  Navier-Stokes when the advection term is implicit and non linearized.
   */

  cs_param_convergence_t        nl_cvg_param;

  /*! \var anderson_param
   * Set of parameters to drive the Anderson acceleration (useful if the type
   * of non-linear algorithm is set to Anderson acceleration).
   */

  cs_iter_algo_param_aac_t      anderson_param;

  /*! @} */

  /*! \var delta_thermal_tolerance
   * Value under which one considers that the thermal equation is converged
   * \f$ max_{c \in Cells} |T_c - T_{ref}|/|T_{ref}| < \epsilon\f$
   * Then stop iteration
   */

  cs_real_t                      delta_thermal_tolerance;

  /*! \var n_max_outer_iter
   * Stopping criterion related to the maximum number of outer iterations
   * allowed. This outer iteration encompasses the Navier-Stokes system,
   * and (according to the case settings) the turbulence system and/or
   * the thermal system.
   */

  int                            n_max_outer_iter;

  /*!
   * @}
   * @name Output
   * @{
   *
   * \var verbosity
   * Level of display of the information related to the Navier-Stokes system
   */

  int                            verbosity;

  /*! \var post_flag
   * Flag storing which predefined post-processing has to be done
   */

  cs_navsto_param_post_flag_t    post_flag;

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

  /*! @}
   * @name Other conditions
   * @{
   */

  /*! \var reference_pressure
   *  Value of the reference pressure p0 (used for rescaling or during update
   *  of physical quantities). By default: 0.
   */

  cs_real_t    reference_pressure;

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
 * \var CS_NSKEY_NL_ALGO
 * Type of algorithm to consider to solve the non-linearity arising from the
 * Navier-Stokes system (Picard or Anderson)
 *
 * \var CS_NSKEY_NL_ALGO_ATOL
 * Absolute tolerance at which the non-linearity arising from the advection
 * term is resolved
 *
 * \var CS_NSKEY_NL_ALGO_DTOL
 * Threshold at which the non-linear algorithm is set as diverged
 *
 * \var CS_NSKEY_NL_ALGO_MAX_ITER
 * Set the maximal number of iterations for solving the non-linearity
 * arising from the advection form
 *
 * \var CS_NSKEY_NL_ALGO_RTOL
 * Relative tolerance at which the non-linearity arising from the advection
 * term is resolved
 *
 * \var CS_NSKEY_SPACE_SCHEME
 * Numerical scheme for the space discretization. Available choices are:
 * - "cdo_fb" or "cdofb" for CDO face-based scheme
 *
 * \var CS_NSKEY_THERMAL_TOLERANCE
 * Value of the tolerance criterion under which one stops iterating on
 * the Navier-Stokes and thermal systems
 *
 * \var CS_NSKEY_VERBOSITY
 * Set the level of details for the specific part related to the Navier-Stokes
 * system
 */

typedef enum {

  CS_NSKEY_DOF_REDUCTION,
  CS_NSKEY_NL_ALGO,
  CS_NSKEY_NL_ALGO_ATOL,
  CS_NSKEY_NL_ALGO_DTOL,
  CS_NSKEY_NL_ALGO_MAX_ITER,
  CS_NSKEY_NL_ALGO_RTOL,
  CS_NSKEY_SPACE_SCHEME,
  CS_NSKEY_THERMAL_TOLERANCE,
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
cs_navsto_param_is_steady(const cs_navsto_param_t       *nsp)
{
  if (nsp == NULL)
    return true;

  if (nsp->model_flag & CS_NAVSTO_MODEL_STEADY)
    return true;
  else
    return false;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure to store all numerical parameters related
 *        to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in] boundaries     pointer to a cs_boundary_t structure
 * \param[in] model          type of model related to the NS system
 * \param[in] model_flag     additional high-level model options
 * \param[in] algo_coupling  algorithm used for solving the NS system
 * \param[in] post_flag      predefined post-processings options
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t          *boundaries,
                       cs_navsto_param_model_t       model,
                       cs_navsto_param_model_flag_t  model_flag,
                       cs_navsto_param_coupling_t    algo_coupling,
                       cs_navsto_param_post_flag_t   post_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a \ref cs_navsto_param_t structure
 *
 * \param[in, out] param  pointer to a \ref cs_navsto_param_t structure
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t *param);

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
 * \brief Apply the numerical settings defined for the Navier-Stokes system to
 *        an equation related to this system. Be aware that the user-defined
 *        settings can be modified in this function when a different choice is
 *        set between the settings for the Navier-Stokes system and the
 *        settings for the momentum equation. The final choice is given by the
 *        settings for the Navier-Stokes system.
 *
 * \param[in]      nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] eqp    pointer to a \ref cs_equation_param_t structure
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
 * \brief Add a new Boussinesq term (source term for the momemtum equation)
 *
 * \param[in, out] nsp               pointer to a cs_navsto_param_t struct.
 * \param[in]      dilatation_coef   value of the dilatation coefficient
 * \param[in]      reference_value   reference value of the associated variable
 *
 * \return a pointer to the newly added structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_boussinesq_t *
cs_navsto_param_add_boussinesq_term(cs_navsto_param_t    *nsp,
                                    cs_real_t             dilatation_coef,
                                    cs_real_t             reference_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the array of values to consider in the Boussinesq term
 *
 * \param[in, out]  bp    pointer to a cs_navsto_param_boussinesq_t structure
 * \param[in]       var   shared pointer to the array of values to consider
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set_boussinesq_array(cs_navsto_param_boussinesq_t   *bp,
                                     const cs_real_t                *var);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         velocity equation (momentum equation in most of the cases)
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the set of parameters related to the momentum equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_navsto_param_get_velocity_param(const cs_navsto_param_t    *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of the model system of equations
 *
 * \param[in]   model    a \ref cs_navsto_param_model_t
 *
 * \return the corresponding name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_model_name(cs_navsto_param_model_t   model);

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
 * \brief Apply the given quadrature rule to all existing definitions under
 *        the cs_navsto_param_t structure
 *
 * \param[in, out] nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in]      qtype    type of quadrature to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set_quadrature_to_all(cs_navsto_param_t    *nsp,
                                      cs_quadrature_type_t  qtype);

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
 * \param[in]      z_name    name of the associated zone (if null or "" if
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
 * \param[in]      z_name    name of the associated zone (if null or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     null or pointer to a structure cast on-the-fly
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
 * \param[in]      z_name    name of the associated zone (if null or "" if
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
 * \param[in]      z_name    name of the associated zone (if null or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     null or pointer to a structure cast on-the-fly
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
 * \brief Add the definition of boundary conditions related to a fixed wall
 *        into the set of parameters for the management of the Navier-Stokes
 *        system of equations
 *
 * \param[in] nsp  pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_fixed_walls(cs_navsto_param_t *nsp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the definition of boundary conditions related to a symmetry
 *        into the set of parameters for the management of the Navier-Stokes
 *        system of equations
 *
 * \param[in] nsp  pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_symmetries(cs_navsto_param_t *nsp);

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
 * \brief  Set the pressure field on a boundary using a uniform value.
 *
 * \param[in] nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in] z_name    name of the associated zone (if null or "" all
 *                      boundary faces are considered)
 * \param[in] values    value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_pressure_bc_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for a sliding wall boundary using a
 *         uniform value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_wall_by_value(cs_navsto_param_t    *nsp,
                                     const char           *z_name,
                                     cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a uniform
 *         value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_value(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_real_t            *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an analytical
 *         function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_analytic(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_analytic_func_t   *ana,
                                         void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an array
 *         of values
 *
 * \param[in]  nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]  z_name       name of the associated zone (if null or "" all
 *                          boundary faces are considered)
 * \param[in]  loc          information to know where are located values
 * \param[in]  array        pointer to an array
 * \param[in]  is_owner     transfer the lifecycle to the cs_xdef_t structure
 *                          (true or false)
 * \param[in]  full_length  if true, the size of "array" should be allocated
 *                          to the total numbers of entities related to the
 *                          given location. If false, a new list is allocated
 *                          and filled with the related subset indirection.
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_array(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_flag_t             loc,
                                      cs_real_t            *array,
                                      bool                  is_owner,
                                      bool                  full_length);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a DoF function
 *
 * \param[in]  nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]  z_name      name of the associated zone (if null or "" all
 *                         boundary faces are considered)
 * \param[in]  dof_loc     where are located DoFs
 * \param[in]  func        pointer to a cs_dof_function_t
 * \param[in]  func_input  null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_dof_func(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_flag_t             dof_loc,
                                         cs_dof_func_t        *func,
                                         void                 *func_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     null or pointer to a structure cast on-the-fly
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
 * \param[in]      z_name    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]      val       value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_val(cs_navsto_param_t *nsp,
                                 const char        *z_name,
                                 cs_real_t          val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an array
 *
 * \param[in] nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in] z_name        name of the associated zone (if null or "" all
 *                          cells are considered)
 * \param[in] loc           information to know where are located values
 * \param[in] array         pointer to an array
 * \param[in] is_owner      transfer the lifecycle to the cs_xdef_t structure
 *                          (true or false)
 * \param[in] full_length   if true, array size is allocated and filled to
 *                          access the full-length array corresponding to
 *                          all locations where are defined the values
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
                                   bool                  full_length);

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
