#ifndef __CS_BOUNDARY_H__
#define __CS_BOUNDARY_H__

/*============================================================================
 * Handle the boundaries of a computational domain
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/* Physic-driven boundary */
typedef enum {

  CS_BOUNDARY_WALL,
  CS_BOUNDARY_SLIDING_WALL,
  CS_BOUNDARY_INLET,
  CS_BOUNDARY_OUTLET,
  CS_BOUNDARY_SYMMETRY,

  /* Physic-driven boundary types for ALE*/
  CS_BOUNDARY_ALE_FIXED,
  CS_BOUNDARY_ALE_SLIDING,
  CS_BOUNDARY_ALE_IMPOSED_VEL,
  CS_BOUNDARY_ALE_IMPOSED_DISP,
  CS_BOUNDARY_ALE_INTERNAL_COUPLING,
  CS_BOUNDARY_ALE_EXTERNAL_COUPLING,
  CS_BOUNDARY_ALE_FREE_SURFACE,

  CS_BOUNDARY_N_TYPES

} cs_boundary_type_t;

/*! \struct cs_boundary_t
 *  \brief Structure storing information related to the "physical" boundaries
 *  that one want to set on the computational domain
 */
typedef struct {

  cs_boundary_type_t    default_type;  /*!< default boundary */

  int                   n_boundaries;  /*!< number of boundaries */
  cs_boundary_type_t   *types;         /*!< type related to each boundary */
  int                  *zone_ids;      /*!< zone id related to each boundary */

} cs_boundary_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_boundary_t  *cs_glob_boundaries; /* Pointer to the shared boundaries
                                            * on the computational domain */

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
cs_boundary_get_name(cs_boundary_type_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to the given \ref cs_boundary_t
 *         structure
 *
 * \param[in, out]   boundaries   pointer to a structure storing boundary info
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_set_default(cs_boundary_t        *boundaries,
                        cs_boundary_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a default boundary structure for the computational domain
 *
 * \param[in]        type         default type of boundary to set
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_boundary_t *
cs_boundary_create(cs_boundary_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 *
 * \param[in, out]   p_boundaries   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_free(cs_boundary_t   **p_boundaries);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in, out] bdy          pointer to a structure storing boundary info
 * \param[in]      type         type of boundary to set
 * \param[in]      zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_add(cs_boundary_t        *bdy,
                cs_boundary_type_t    type,
                const char           *zone_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_BOUNDARY_WALL zone type
 *
 * \param[in, out]  boundaries    pointer to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_def_wall_zones(cs_boundary_t   *boundaries);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 *
 * \param[in] bdy          pointer to a structure storing boundary info
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_log_setup(const cs_boundary_t     *bdy);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_H__ */
