#ifndef __CS_MESH_DEFORM_H__
#define __CS_MESH_DEFORM_H__

/*============================================================================
 * Mesh deformation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 * \param[in]        n_zones    number of zones with prescribed deformation
 * \param[in]        zones_ids  ids of zones with prescribed deformation
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_initialize(cs_domain_t   *domain,
                          int            n_zones,
                          const int      zone_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribe the displacement vector for a set of vertices.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used.
 *
 * This function may be called multiple times before displacements are
 * are computed, so displacements may be prescribed for different zones
 * separately or simultaneously.
 *
 * \param[in]  n_vertices    number of vertices at which to prescribe
 *                           displacements
 * \param[in]  vertex_ids    ids of vertices at which to prescribe
 *                           displacements, or NULL if [0, ... n_vertices-1]
 * \param[in]  displacement  pointer to prescribed displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_prescribe_displacement(cs_lnum_t          n_vertices,
                                      const cs_lnum_t    vertex_ids[],
                                      const cs_real_3_t  displacement[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute displacement for mesh deformation.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_solve_displacement(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to current mesh displacement vector.
 *
 * \return  pointer to current displacement vector
 */
/*----------------------------------------------------------------------------*/

const cs_real_3_t *
cs_mesh_deform_get_displacement(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free structures used fo mesh deformation.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_DEFORM_H__ */
