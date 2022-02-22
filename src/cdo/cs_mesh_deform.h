#ifndef __CS_MESH_DEFORM_H__
#define __CS_MESH_DEFORM_H__

/*============================================================================
 * Mesh deformation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_domain.h"
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
 * \brief  Test if mesh deformation is activated
 *
 * \return true if mesh deformation computation is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_mesh_deform_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future mesh deformation
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_activate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the boundary zones on which mesh deformation is prescribed.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used.
 *
 * \param[in]  n_boundary_zones   number of boundary zones at which to
 *                                prescribe displacements
 * \param[in]  boundary_zone_ids  ids of boundary zones at which to
 *                                prescribe displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_define_dirichlet_bc_zones(cs_lnum_t  n_boundary_zones,
                                         const int  boundary_zone_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup the equations related to mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_setup(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribe the displacement vector for a set of vertices.
 *
 * Only those values at vertices matching boundary zones with prescribed
 * displacement will really be used, as defined by
 * \ref cs_mesh_deform_define_dirichlet_bc_zones.
 *
 * When calling this function multiple times for different vertex sets,
 * the most recently defined values are used for vertices belonging to
 * multiple sets.
 *
 * \param[in]  n_vertices         number of vertices at which to prescribe
 *                                displacements
 * \param[in]  vertex_ids         ids of vertices at which to prescribe
 *                                displacements, or NULL if
 *                                [0, ... n_vertices-1]
 * \param[in]  displacement       pointer to prescribed displacements
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_prescribe_displacement(cs_lnum_t          n_vertices,
                                      const cs_lnum_t    vertex_ids[],
                                      const cs_real_3_t  displacement[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a fixed displacement vector for given vertices.
 *
 * This displacement is enforced at all given vertices, including interior
 * vertices.
 *
 * If this function is called multiple times, the previous definitions
 * are overwritten, so all displacements of this type must be defined
 * in a single call to this function.
 *
 * \param[in]  n_vertices         number of vertices at which to prescribe
 *                                displacements
 * \param[in]  vertex_ids         ids of vertices at which to prescribe
 *                                displacements, or NULL if
 *                                [0, ... n_vertices-1]
 * \param[in]  displacement       pointer to prescribed displacements,
 *                                or NULL for no displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_force_displacements(cs_lnum_t          n_vertices,
                                   const cs_lnum_t    vertex_ids[],
                                   const cs_real_3_t  displacement[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute displacement for mesh deformation.
 *
 * \param[in, out]   domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_deform_solve_displacement(cs_domain_t  *domain);

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
