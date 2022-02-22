#ifndef __CS_MESH_BOUNDARY_H__
#define __CS_MESH_BOUNDARY_H__

/*============================================================================
 * Insert boundaries into the mesh.
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
 * \brief Insert boundary into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 * The created faces share vertices.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \param[in, out]  mesh     pointer to mesh structure to modify
 * \param[in]       n_faces  number of selected (interior) faces
 * \param[in, out]  face_id  list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert(cs_mesh_t  *mesh,
                        cs_lnum_t   n_faces,
                        cs_lnum_t   face_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundary into the mesh, sharing vertices on both sides.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * This function should be used with caution, as sharing vertices
 * may cause issues with vertex-based values, gradient extended
 * neighborhoods, and some visualization operations.
 *
 * \param[in, out]  mesh     pointer to mesh structure to modify
 * \param[in]       n_faces  number of selected (interior) faces
 * \param[in, out]  face_id  list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert_with_shared_vertices(cs_mesh_t  *mesh,
                                             cs_lnum_t   n_faces,
                                             cs_lnum_t   face_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert a boundary into the mesh between a given set of cells
 *        and the the others.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * An optional group name may be assigned to newly created boundary faces.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * \param[in, out]  mesh         pointer to mesh structure to modify
 * \param[in]       group_name  group name to assign to newly created boundary
 *                              faces, or NULL
 * \param[in]       n_cells     number of separated cells
 * \param[in]       cell_id     separated cell ids
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert_separating_cells(cs_mesh_t        *mesh,
                                         const char       *group_name,
                                         cs_lnum_t         n_cells,
                                         const cs_lnum_t   cell_id[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove periodicity information from a mesh.
 *
 * Periodic interior faces are transformed into boundary faces.
 *
 * \param[in, out]  mesh  pointer to mesh structure to modify
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_remove_periodicity(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BOUNDARY_H__ */
