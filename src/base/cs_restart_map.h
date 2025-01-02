#ifndef __CS_RESTART_MAP_H__
#define __CS_RESTART_MAP_H__

/*============================================================================
 * Checkpoint / restart extension to mapped meshes
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Local type definitions
 *============================================================================*/


/*=============================================================================
 * Global variables
 *============================================================================*/


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate restart files should be mapped to a given mesh input
 *
 * \param[in]  mesh_path           path to mesh input
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_mesh_input(const char  *mesh_path);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options relative to restart file mapping to a given mesh input.
 *
 * \param[in]  apply_mesh_deformation  apply mesh deformation from upstream
 *                                     computation (if present) so as to map
 *                                     to final, and not initial mesh shape.
 * \param[in]  tolerance_base          associated base tolerance (used for
 *                                     bounding box check only, not for
 *                                     location test)
 * \param[in]  tolerance_fraction      associated fraction of element bounding
 *                                     boxes added to tolerance
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_options(bool   apply_mesh_deformation,
                           float  tolerance_base,
                           float  tolerance_fraction);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate whether location for restart file mapping is needed at
 *         cells or vertices.
 *
 * By default, mapping is done for cell-based quantities, but not for
 * vertex-based quantities.
 *
 * Mapping of quantities at faces or particles is not handled yet, but will
 * use the cell-center or vertex based mappings in the future in all cases:
 * - interior faces may not be aligned with previous faces, so some sort of
 *   cell-based interpolation will be required
 * - boundary faces can use the boundary face / cell adjacency to avoid an
 *   additional mapping
 * - for particles, as the previous location is stored based on cell ids,
 *   updating particle locations will require locating them in the matching
 *   cell and completing a trajectory adjustment (direct location should be
 *   avoided in case of concave boundaries).
 *
 * \param[in]  map_cell_centers    locate cell centers in the previous mesh.
 * \param[in]  map_vertices        locate vertices in the previous mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_set_locations(bool map_cell_centers,
                             bool map_vertices);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build mapping of restart files to different mesh if defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_build(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free restart file mapping to different mesh.
 *
 * Revert restart reading to default behavior.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_map_free(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_MAP_H__ */
