#ifndef __CS_NAVSTO_PARAM_H__
#define __CS_NAVSTO_PARAM_H__

/*============================================================================
 * Routines to handle cs_navsto_param_t structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_xdef.h"

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
  CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES,
  CS_NAVSTO_MODEL_BOUSSINESQ_NAVIER_STOKES,

  CS_NAVSTO_N_MODELS

} cs_navsto_param_model_t;

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
 * \var CS_NAVSTO_COUPLING_UZAWA
 * The system is solved without decoupling the equations using a Uzawa algorithm
 * and an Augmented Lagrangian approach inside each sub-iteration.
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
 * \var CS_NAVSTO_COUPLING_PROJECTION
 * The system is solved using an incremental projection algorithm
 */

typedef enum {

  CS_NAVSTO_COUPLING_UZAWA,
  CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY,
  CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP,
  CS_NAVSTO_COUPLING_PROJECTION,

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
  int                           verbosity;

  /*! \var dof_reduction_mode
   *  How are defined the Degrees of freedom
   */
  cs_param_dof_reduction_t      dof_reduction_mode;

  /*! \var time_scheme
   * Discretization scheme for time
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
 * \var CS_NSKEY_QUADRATURE
 * Set the type to use in all routines involving quadrature (similar to \ref
 * CS_EQKEY_BC_QUADRATURE)
 *
 * \var CS_NSKEY_SPACE_SCHEME
 * Numerical scheme for the space discretization
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
  CS_NSKEY_QUADRATURE,
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
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in]  model          model related to the NS system to solve
 * \param[in]  time_state     state of the time for the NS equations
 * \param[in]  algo_coupling  algorithm used for solving the NS system
*
 * \return a pointer to a new allocated structure
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
 * \param[in]  model          model related to the NS system to solve
 * \param[in]  time_state     state of the time for the NS equations
 * \param[in]  algo_coupling  algorithm used for solving the NS system
*
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(cs_navsto_param_model_t        model,
                       cs_navsto_param_time_state_t   time_state,
                       cs_navsto_param_coupling_t     algo_coupling);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a cs_navsto_param_t
 *         structure
 *
 * \param[in, out] nsp      pointer to a cs_navsto_param_t structure to set
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
                                   cs_lnum_t            *index);

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
 * \param[in]     coupling    A \ref cs_navsto_param_coupling_t
 *
 * \return the name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_coupling_name(cs_navsto_param_coupling_t  coupling);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_PARAM_H__ */
