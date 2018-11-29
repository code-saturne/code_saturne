#ifndef __CS_BOUNDARY_H__
#define __CS_BOUNDARY_H__

/*============================================================================
 * Handle the boundaries of a computational domain
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
#define CS_BOUNDARY_WALLS_NAME   "cs_boundary_walls"

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Physic-driven boundary types */
typedef enum {

  CS_BOUNDARY_WALL,
  CS_BOUNDARY_SLIDING_WALL,
  CS_BOUNDARY_INLET,
  CS_BOUNDARY_OUTLET,
  CS_BOUNDARY_SYMMETRY,
  CS_DOMAIN_N_BOUNDARY_TYPES

} cs_boundary_type_t;

/* Structure for handling the boundaries of the domain */
typedef struct {

  cs_boundary_type_t    default_type;

  int                   n_boundaries;
  cs_boundary_type_t   *types;
  int                  *zone_ids;

} cs_boundary_t;

extern cs_boundary_t  *cs_glob_boundaries;

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
cs_boundaries_get_name(cs_boundaries_type_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to this domain
 *
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_boundaries_set_default(cs_boundaries_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the default boundary associated to the computational domain
 *
 * \return the type of the domain boundary defined by default
 */
/*----------------------------------------------------------------------------*/

cs_boundaries_type_t
cs_boundaries_get_default(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_boundaries_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in] type         type of boundary to set
 * \param[in] zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_boundaries_add(cs_boundaries_type_t    type,
                       const char                  *zone_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_BOUNDARIES_WALL zone type
 */
/*----------------------------------------------------------------------------*/

void
cs_boundaries_def_wall_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_boundaries_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARIES_H__ */
