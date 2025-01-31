#ifndef __CS_BOUNDARY_CONDITIONS_PRIV_H__
#define __CS_BOUNDARY_CONDITIONS_PRIV_H__

/*============================================================================
 * Private types for boundary condition handling.
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

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_function.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_time_control.h"
#include "base/cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*!< Velocity normalization modes */

typedef enum {

  CS_BC_VEL_RESCALE_NONE,             /*!< value is direct velocity */
  CS_BC_VEL_RESCALE_MASS_FLOW_RATE,   /*!< value represents mass flow rate */
  CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE  /*!< value represents volume flow rate */

} cs_bc_velocity_rescale_t;

/*!< Turbulence computation modes */

/* Enum for boundary conditions of multiple turbulence values */

typedef enum {

  CS_BC_TURB_NONE = 0,                /*!< no computation with this mechanism */
  CS_BC_TURB_BY_HYDRAULIC_DIAMETER,   /*!< based on hydraulic diameter */
  CS_BC_TURB_BY_TURBULENT_INTENSITY   /*!< based on turbulent intensity */

} cs_bc_turbulence_compute_t;

/*! the given value is constant in time */
#define CS_BC_OPEN_CONSTANT           (1 << 1)

/*! the given quantity is uniform over the zone */
#define CS_BC_OPEN_UNIFORM_QUANTITY   (1 << 2)

/*! the velocity vector is normal to the boundary */
#define CS_BC_OPEN_UNIFORM_DIRECTION  (1 << 3)

/*! the velocity vector is normal to the boundary */
#define CS_BC_OPEN_NORMAL_DIRECTION   (1 << 4)

/*! \brief Inlet definition (velocity and turbulence) */
/*----------------------------------------------------*/

/*! Since inlets (and outlets) may involve coupling several variables,
  a common context may be needed for xdef-based BCs for those variables.

  Also, in cases where MPI or thread synchronization may be needed,
  xdef-based definitions might not be the most appropriate, or need
  "by array/by field" type semantics, (or use "dof" function with
  pre-synchronized data).

  The following structure allows associating 2 steps to the inlet
  velocity computation:

  - Computing the local velocity or its norm.
  - Optionally, apply a scaling function.

  Performance may be slightly better if scaling is done direclty in the first
  function (if relevant), but allowing an optional second step allows
  better mapping to some specific models.
*/

typedef struct {

  const  cs_zone_t       *zone;               /*!< Pointer to zone */

  cs_time_control_t       tc;                 /*!< Time control sub-structure */

  int                     vel_flags;          /*!< input mode flags */

  cs_bc_velocity_rescale_t    vel_rescale;    /*!< velocity rescaling type */
  cs_bc_turbulence_compute_t  turb_compute;   /*!< turbulence computation mode */

  int                     bc_pm_zone_num;     /*!< matching zone number
                                                in cs_glob_bc_pm_info;
                                                (legacy mapping) */

  cs_real_t               vel_values[4];      /*!< velocity vector and value or
                                                flow rate when uniform
                                                [u_x, u_y, u_z, 0] or
                                                [0, 0, 0, u_norm] or
                                                [x, y, z, q] */

  cs_real_t              *vel_buffer;         /*!< Pointer to per face or uniform
                                                velocity (depending on method) */

  cs_real_t               hyd_diameter;       /*!< hydraulic diameter if > 0 */
  cs_real_t               turb_intensity;     /*!< turbulence intensity if > 0 */

  cs_real_t               c_pr;               /*!< imposed pressure
                                                (compressible flow only) */
  cs_real_t               c_tk;               /*!< imposed temperature
                                                (compressible flow only) */

  cs_eval_at_location_t  *vel_func;           /*!< associated vector evaluation
                                                function, or NULL */
  void                   *vel_func_input;     /*!< Optional vector evaluation
                                                input, or NULL */

  cs_eval_at_location_t  *flow_func;          /*!< associated global (mass or
                                                volume) flow evaluation
                                                function, or NULL */
  void                   *flow_func_input;    /*!< Optional flow evaluation
                                                evaluation input, or NULL */

  cs_eval_at_location_t  *scale_func;         /*!< Associated scaling
                                                function, or NULL */
  void                   *scale_func_input;   /*!< Optional scaling function
                                                input, or NULL */

  cs_dof_func_t          *dof_func;           /*!< Associated dof_func
                                                (keep track of it for simpler
                                                tests) */

  void                   *model_inlet;        /*!< Additional model-based
                                                inlet structure */

  cs_destructor_t        *model_inlet_del;    /*!< Destructor for associated
                                                inlet structure (use simple
                                                CS_FREE if no destructor
                                                is associated) */

} cs_boundary_conditions_open_t;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find or add an open boundary context for a given zone.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_boundary_conditions_open_t *
cs_boundary_conditions_open_find_or_add(const  cs_zone_t   *zone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an open boundary context structure for a given zone.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone, or NULL if not found.
 */
/*----------------------------------------------------------------------------*/

cs_boundary_conditions_open_t *
cs_boundary_conditions_open_find(const  cs_zone_t   *zone);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_PRIV_H__ */
