#ifndef __CS_CDOFB_MONOLITHIC_PRIV_H__
#define __CS_CDOFB_MONOLITHIC_PRIV_H__

/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
 * coupling)
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_monolithic_sles.h"
#include "cs_cdofb_navsto.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_navsto_coupling.h"
#include "cs_param.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic_priv.c
 *
 * \brief Structures and function pointers useful to build and solve
 *        the Navier-Stokes equations with face-based schemes and a
 *        monolithic approach
 */

typedef struct _cdofb_monolithic_t  cs_cdofb_monolithic_t;

/*=============================================================================
 * Definitions of function pointers
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *
 * \param[in]       csys              pointer to a cs_cell_sys_t structure
 * \param[in]       cm                pointer to a cs_cell_mesh_t structure
 * \param[in]       div_op            array with the divergence op. values
 * \param[in]       has_sourceterm    has the equation a source term?
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  eqa               pointer to cs_equation_assemble_t
 * \param[in, out]  mav               pointer to cs_matrix_assembler_values_t
 * \param[in, out]  rhs               right-end side of the system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_monolithic_assemble_t)(const cs_cell_sys_t            *csys,
                                 const cs_cell_mesh_t           *cm,
                                 const cs_real_t                *div_op,
                                 const bool                      has_sourceterm,
                                 cs_cdofb_monolithic_t          *sc,
                                 cs_cdofb_vecteq_t              *eqc,
                                 cs_equation_assemble_t         *eqa,
                                 cs_matrix_assembler_values_t   *mav,
                                 cs_real_t                       rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         steady-state case. Specific case: GKB algorithm is used to solve
 *         the saddle-point system.
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in]      dir_values    array storing the Dirichlet values
 * \param[in]      forced_ids    indirection in case of internal enforcement
 * \param[in, out] sc            pointer to the scheme context
 * \param[in, out] matrix        pointer to a \ref cs_matrix_t structure
 * \param[in, out] mom_rhs       rhs array related to the momentum eq.
 * \param[in, out] mass_rhs      rhs array related to the mass eq.
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_monolithic_build_t)(const cs_navsto_param_t      *nsp,
                              const cs_real_t              *dir_values,
                              const cs_lnum_t               forced_ids[],
                              cs_cdofb_monolithic_t        *sc,
                              cs_matrix_t                  *matrix,
                              cs_real_t                    *mom_rhs,
                              cs_real_t                    *mass_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

typedef int
(cs_cdofb_monolithic_solve_t)(const cs_navsto_param_t       *nsp,
                              const cs_equation_param_t     *eqp,
                              cs_cdofb_monolithic_sles_t    *msles);

/*=============================================================================
 * Structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_monolithic_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         Navier-Stokes equations and vector-valued face unknowns.
 *         Case of a monolithic approach (i.e fully coupled)
 */

struct _cdofb_monolithic_t {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_monolithic_t (owned by \ref
   *  cs_navsto_system_t) containing the settings related to the monolithic
   *  approach
   */

  cs_navsto_monolithic_t   *coupling_context;

  /*!
   * @name Main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   */

  /*! \var velocity
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the velocity
   */

  cs_field_t               *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t               *pressure;

  /*! \var divergence
   *  Pointer to \ref cs_real_t containing the values of the divergence on the
   *  cells
   */

  cs_field_t               *divergence;

  /*!
   * @}
   * @name Boundary conditions (BC) management
   * Routines and elements used for enforcing the BCs
   * @{
   *
   *  \var bf_type
   *  Array of boundary type for each boundary face. (Shared)
   */

  const cs_boundary_type_t       *bf_type;

  /*!
   * \var pressure_bc
   * Structure storing the metadata after processing the user-defined boundary
   * conditions related to the pressure field
   */

  cs_cdo_bc_face_t               *pressure_bc;

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

  cs_cdo_apply_boundary_t        *apply_fixed_wall;
  cs_cdo_apply_boundary_t        *apply_sliding_wall;
  cs_cdo_apply_boundary_t        *apply_velocity_inlet;
  cs_cdo_apply_boundary_t        *apply_symmetry;

  /*!
   * @}
   * @name Build stage
   * Additional members which corresponds to function pointers
   * @{
   */

  cs_cdofb_monolithic_build_t      *steady_build;
  cs_cdofb_monolithic_build_t      *build;
  cs_cdofb_monolithic_assemble_t   *assemble;

  /*!
   * @}
   * @name Solve stage
   * Additional members which may be used to solve the system
   * @{
   */

  cs_cdofb_monolithic_solve_t      *solve;

  /* \var msles
   * Set of pointers to enable the resolution of saddle-point system
   * with various algorithms. This structure allows us to unify the prototype
   * of "solve" functions
   * Some members of this structure are allocated only if a specific algorithm
   * is requested.
   */

  cs_cdofb_monolithic_sles_t       *msles;

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
  cs_timer_counter_t  timer;

  /*! @} */

};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_MONOLITHIC_PRIV_H__ */
