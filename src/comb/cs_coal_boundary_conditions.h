#ifndef __CS_COAL_BOUNDARY_CONDITIONS_H__
#define __CS_COAL_BOUNDARY_CONDITIONS_H__

/*============================================================================
 * Coal combustion model boundary conditions.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_coal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*! \brief Inlet model state data */

typedef struct _cs_coal_bc_inlet_state_t  cs_coal_bc_inlet_state_t;

/*! \brief Inlet definition for pulverized coal combustion */

typedef struct {

  const  cs_zone_t  *zone;                     /*!< Pointer to associated zone */

  int  ientat;                                 /*!< Air indicator */
  int  ientcp;                                 /*!< Cp indicator */
  int  inmoxy;                                 /*!< Number of oxydants */

  cs_real_t  qimpcp[CS_COMBUSTION_MAX_COALS];  /*!< Coal mass flow per coal */
  cs_real_t  timpcp[CS_COMBUSTION_MAX_COALS];  /*!< Coal temperature (in K)
                                                 per coal */

  /*! Coal class mass distribution ratio (in percentage) */

  cs_real_t  distch[CS_COMBUSTION_MAX_COALS]
                   [CS_COMBUSTION_MAX_CLASSES_PER_COAL];

  cs_real_t  t_air;                            /*!< Air temperature (in K) */

  /*! @name Prescribed air mass flow definition.
   * @{

     Using such a definition allows prescribing
    the air mass flow and letting the code determine the total
    (air + coal) mass flow.

    Only used if activated by call to
    \ref cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_value or
    \ref cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_func.
  */

  cs_real_t               qm_air;        /*!< Air mass flow rate,
                                           if prescribed */
  cs_eval_at_location_t  *qm_air_func;   /*!< Associated global air mass flow
                                           evaluation function, or NULL */
  void                   *qm_air_input;  /*!< Associated function
                                           input if needed, or NULL */

  /*! @} */

  cs_coal_bc_inlet_state_t  *state;      /* Inlet state (not user-accessible) */

} cs_coal_bc_inlet_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to coal boundary conditions inlet structure.
 *
 * If no such structure was previously present, it is created and linked
 * to the matching open boundary condition inlet.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_coal_bc_inlet_t *
cs_coal_boundary_conditions_get_inlet(const  cs_zone_t   *zone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant air mass flow rate to an inlet.
 *
 * The total mass flow rate will also include that of the pulverized coals.
 *
 * This is otherwise similar to
 * \ref cs_boundary_conditions_open_set_mass_flow_rate_by_value.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_value
  (const  cs_zone_t  *z,
   cs_real_t          q);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign an air mass flow rate to an inlet based on provided function.
 *
 * The total mass flow rate will also include that of the pulverized coals.
 *
 * This is otherwise similar to
 * \ref cs_boundary_conditions_open_set_mass_flow_rate_by_func.
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (mass flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for pulverized coal combution.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density at inlets for pulverized coal combustion.
 *
 * This is based on boundary condition definitions, but is called at an
 * earlier stage in the time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_boundary_conditions_inlet_density(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COAL_BOUNDARY_CONDITIONS_H_ */
