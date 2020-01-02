#ifndef __FVM_NODAL_EXTRUDE_H__
#define __FVM_NODAL_EXTRUDE_H__

/*============================================================================
 * Extrusion of a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extrude nodal mesh.
 *
 * Vertex and element parent numbering is removed if present.

 * Note: layout of new elements in memory is such that the definitions
 *       of all elements extruded from a same ancestor are contiguous.
 *       that is, {e_1, e_2, ..., e_n} leads to
 *       {e_1_layer_1, ..., e_1_layer_m, e_2_layer_1, ... e_n_layer_m}
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be extruded
 *   n_layers          <-> number of extruded layers
 *   extrusion_vectors <-> length and direction of extrusion for each vertex;
 *                         size: mesh_spatial_dim . n_vertices
 *   distribution      <-> optional distribution of resulting vertices
 *                         along each extrusion vector (size: n_layers + 1)
 *                         with values ranging from 0 to 1, or NULL for
 *                         a regular distribution.
 *----------------------------------------------------------------------------*/

void
fvm_nodal_extrude(fvm_nodal_t        *this_nodal,
                  const cs_lnum_t     n_layers,
                  const cs_coord_t    extrusion_vectors[],
                  const cs_coord_t    distribution[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_EXTRUDE_H__ */
