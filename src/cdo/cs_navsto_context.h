#ifndef CS_NAVSTO_CONTEXT_H
#define CS_NAVSTO_CONTEXT_H

/*============================================================================
 * Functions to handle cs_navsto_system_t structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_field.h"
#include "cdo/cs_equation_bc.h"
#include "cdo/cs_iter_algo.h"
#include "cdo/cs_navsto_coupling.h"

namespace cs {

/*=============================================================================
 * Structure definitions
 *============================================================================*/

/*! \struct cdo_navsto_ctx_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         Navier-Stokes equations and vector-valued face unknowns.
 */

struct cdo_navsto_ctx_t {
  /*!
   * @name Main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   */

  /*! \var velocity
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the velocity
   */

  cs_field_t *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t *pressure;

  /*! \var divergence
   *  Pointer to \ref cs_real_t containing the values of the divergence on the
   *  cells
   */

  cs_field_t *divergence;

  /*!
   * @}
   * @name Advection quantities
   * Members related to the advection
   * @{
   *
   *  \var adv_field
   *  Pointer to the cs_adv_field_t related to the Navier-Stokes eqs (Shared)
   */

  cs_adv_field_t *adv_field;

  /*! \var mass_flux_array
   *  Current values of the mass flux at primal faces (Shared)
   */

  cs_real_t *mass_flux_array;

  /*! \var mass_flux_array_pre
   *  Previous values of the mass flux at primal faces (Shared)
   */

  cs_real_t *mass_flux_array_pre;

  /*!
   * @}
   * @name Boundary conditions (BC) management
   * Functions and elements used for enforcing the BCs
   * @{
   *
   *
   *  \var bf_type
   *  Array of boundary type for each boundary face. (Shared)
   */

  const cs_boundary_type_t *bf_type;

  /*!
   * \var pressure_bc
   * Structure storing the metadata after processing the user-defined boundary
   * conditions related to the pressure field
   */

  cs_cdo_bc_face_t *pressure_bc;
  int               pressure_rescaling;

  /*! \var apply_fixed_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (no slip boundary)
   *
   *  \var apply_sliding_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (a tangential velocity is specified at the wall)
   *
   *  \var apply_velocity_inlet
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  boundary with a fixed velocity at the inlet
   *
   *  \var apply_symmetry
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  symmetry boundary
   */

  cs_cdo_apply_boundary_t *apply_fixed_wall;
  cs_cdo_apply_boundary_t *apply_sliding_wall;
  cs_cdo_apply_boundary_t *apply_velocity_inlet;
  cs_cdo_apply_boundary_t *apply_symmetry;

  /*!
   * @}
   * @name Performance monitoring
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   * @{
   */

  /*! \var timer
   *  Cumulated elapsed time for building and solving the Navier--Stokes system
   */

  cs_timer_counter_t timer;

  /*! @} */

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Retrieve the mass flux array related to the Navier-Stokes system.
   *
   * \param[in] previous    if true return the previous state otherwise the
   *                        current state.
   *
   * \return a pointer to an array of cs_real_t
   */
  /*----------------------------------------------------------------------------*/

  cs_real_t *
  get_mass_flux(bool previous) const
  {
    cs_real_t *mass_flux = mass_flux_array;
    if (previous)
      mass_flux = mass_flux_array_pre;

    return mass_flux;
  }
};

/*! \struct cdo_navsto_monolithic_ctx_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         Navier-Stokes equations and vector-valued face unknowns.
 *         Case of a monolithic approach (i.e fully coupled)
 */

struct cdo_navsto_monolithic_ctx_t : public cdo_navsto_ctx_t {
  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_monolithic_t (owned by \ref
   *  cs_navsto_system_t) containing the settings related to the monolithic
   *  approach
   */

  cs_navsto_monolithic_t *coupling_context;
};

/*! \struct cdo_navsto_predco_ctx_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         Navier-Stokes equations with a prediction/correction algorithm
 */

struct cdo_navsto_predco_ctx_t : public cdo_navsto_ctx_t {
  /*! \var coupling_context

   *  Pointer to a \ref cs_navsto_projection_t_t (owned by
   *  \ref cs_navsto_system_t) containing the settings related to a prjection
   *  or prediction/correction algorithm.
   */

  cs_navsto_projection_t *coupling_context;

  /*!
   * @name Additional main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   * @{
   */

  /*! \var predicted_velocity_f
   * Values of the predicted velocity at faces.
   * This values may not be divergence-free
   */
  cs_real_t *predicted_velocity_f;

  /*! @} */
};

/*! \var coupling_context
 *  Pointer to a \ref cdo_navsto_ac_ctx_t (owned by \ref cs_navsto_system_t)
 *  containing the settings related to an artificial compressibility (AC)
 *  algorithm
 */

struct cdo_navsto_ac_ctx_t : public cdo_navsto_ctx_t {
  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_ac_t (owned by \ref cs_navsto_system_t)
   *  containing the settings related to an artificial compressibility (AC)
   *  algorithm
   */

  cs_navsto_ac_t *coupling_context;

  /*!
   *
   * @name Parameters of the algorithm
   * Easy access to useful features and parameters of the algorithm
   * @{
   */

  /*! \var is_zeta_uniform
   *  Bool telling if the auxiliary parameter zeta is uniform. Not always
   *  necessary: zeta is typically used in Artificial Compressibility algos
   */

  bool is_zeta_uniform;

  /*!
   * @}
   * @name Convergence monitoring
   * Structure used to drive the convergence of high-level iterative algorithms
   * @{
   */

  /*!
   * \var nl_algo
   * Structure driving the convergence of the non-linear algorithm
   */

  cs_iter_algo_t *nl_algo;

  /*! @} */
};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

} // end namespace cs

#endif /* CS_NAVSTO_CONTEXT_H */
