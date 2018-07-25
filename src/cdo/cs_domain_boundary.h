#ifndef __CS_DOMAIN_BOUNDARY_H__
#define __CS_DOMAIN_BOUNDARY_H__

/*============================================================================
 * Handle the boundary of a computational domain
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

#include "cs_base.h"
#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Name of the boundary zone gathering all domain boundary walls */
#define CS_DOMAIN_BOUNDARY_WALLS_NAME   "cs_domain_boundary_walls"

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Physic-driven boundary */
typedef enum {

  CS_DOMAIN_BOUNDARY_WALL,
  CS_DOMAIN_BOUNDARY_SLIDING_WALL,
  CS_DOMAIN_BOUNDARY_INLET,
  CS_DOMAIN_BOUNDARY_OUTLET,
  CS_DOMAIN_BOUNDARY_SYMMETRY,
  CS_DOMAIN_N_BOUNDARY_TYPES

} cs_domain_boundary_type_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the domain boundary condition
 *
 * \param[in] type     type of domain boundary
 *
 * \return the associated boundary name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_domain_boundary_get_name(cs_domain_boundary_type_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to this domain
 *
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_set_default(cs_domain_boundary_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the default boundary associated to the computational domain
 *
 * \return the type of the domain boundary defined by default
 */
/*----------------------------------------------------------------------------*/

cs_domain_boundary_type_t
cs_domain_boundary_get_default(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in] type         type of boundary to set
 * \param[in] zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_add(cs_domain_boundary_type_t    type,
                       const char                  *zone_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_DOMAIN_BOUNDARY_WALL zone type
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_def_wall_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_boundary_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_BOUNDARY_H__ */
