#ifndef __FVM_INTERPOLATE_H__
#define __FVM_INTERPOLATE_H__

/*============================================================================
 * Interpolate data defined on a nodal representation associated with a mesh.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm/fvm_defs.h"
#include "fvm/fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Interpolate vertex-based values to points located relative to a mesh.
 *
 * Least-squared based interpolation is used for now.
 *
 * parameters:
 *   this_nodal         <-- pointer to nodal mesh representation structure
 *   entity_dim         <-- 3 for location on cells, 2 for faces, 1 for edges,
 *                          and 0 for vertices
 *   data_dim           <-- data dimension
 *   n_points           <-- number of points to locate
 *   location_id        <-- id of element (with concatenated sections)
 *                          in which each point is located
 *   point_coords       <-- point coordinates
 *   src_data           <-- source data (interleaved)
 *   dest_data          <-> destination data (interleaved)
 *----------------------------------------------------------------------------*/

void
fvm_interpolate_vtx_data(const fvm_nodal_t       *this_nodal,
                         int                      entity_dim,
                         int                      data_dim,
                         cs_lnum_t                n_points,
                         const cs_lnum_t          location_id[],
                         const cs_coord_t         point_coords[],
                         const cs_real_t          src_data[],
                         cs_real_t                dest_data[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_INTERPOLATE_H__ */
