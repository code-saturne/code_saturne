#ifndef CS_MESH_CUT_H
#define CS_MESH_CUT_H

/*============================================================================
 * Functions to remove mesh elements.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "mesh/cs_mesh.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Cut cells with planes.
 *
 * Each cell can be cut by a single plane, defined by its normal and an origin
 * (i.e. any point in the plane). Cells whose assigned normals are null
 * vectors are not cut.
 *
 * The polygons created by the cut are added to a new group,
 * "auto:closing_polygons".
 *
 * This function should be followed by applying a joining on the group,
 * "auto:transformed_internal_faces".
 *
 * \param[in, out]  mesh      mesh to cut
 * \param[in]       p_normals array of plane_normals of size mesh->n_cells
 * \param[in]       p_origins array of plane origins of size mesh->n_cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cut(cs_mesh_t       *mesh,
            const cs_real_t  p_normals[][3],
            const cs_real_t  p_origins[][3]);

/*----------------------------------------------------------------------------*/

#endif /* CS_MESH_CUT_H */
