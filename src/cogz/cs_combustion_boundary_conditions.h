#ifndef CS_COMBUSTION_BOUNDARY_CONDITIONS_H
#define CS_COMBUSTION_BOUNDARY_CONDITIONS_H

/*============================================================================
 * Gas combustion model boundary conditions.
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

#include "base/cs_defs.h"
#include "cogz/cs_combustion_gas.h"
#include "base/cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*! \brief Inlet definition for pulverized combustion combustion */

typedef struct {

  const  cs_zone_t  *zone;  /*!< Pointer to associated zone */

  int  ientox;              /*!< Oxydant inlet indicator */
  int  ientfu;              /*!< Fuel inlet indicator */
  int  ientgf;              /*!< Fresh gas indicator */
  int  ientgb;              /*!< Burned gas indicator */

  cs_real_t  fment;         /*!< Mixture fraction */
  cs_real_t  tkent;         /*!< Inlet temperaturen in K */

  cs_real_t  tgf;           /*!< Inlet temperaturen in K */

} cs_combustion_bc_inlet_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to combustion boundary conditions inlet structure.
 *
 * If no such structure was previously present, it is created and linked
 * to the matching open boundary condition inlet.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_combustion_bc_inlet_t *
cs_combustion_boundary_conditions_get_inlet(const  cs_zone_t   *zone);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion with EBU model.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_ebu(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion with Libby-Williams
 *        model.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_lw(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density at boundary for pulverized combustion combustion.
 *
 * This is based on boundary condition definitions, but is called at an
 * earlier stage in the time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_density(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COMBUSTION_BOUNDARY_CONDITIONS_H */
