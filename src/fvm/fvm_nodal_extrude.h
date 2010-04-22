#ifndef __FVM_NODAL_EXTRUDE_H__
#define __FVM_NODAL_EXTRUDE_H__

/*============================================================================
 * Extrusion of a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2005-2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
 *   this_section      <-> pointer to structure that should be extruded
 *   n_layers          <-> number of extruded layers
 *   extrusion_vectors <-> length and direction of extrusion for each vertex;
 *                         size: mesh_spatial_dim . n_vertices
 *   distribution      <-> optional distribution of resulting vertices
 *                         along each extrusion vector (size: n_layers + 1)
 *                         with values ranging from 0 to 1, or NULL for
 *                         a regular distribution.
 *----------------------------------------------------------------------------*/

void
fvm_nodal_extrude(fvm_nodal_t        *const this_nodal,
                  const fvm_lnum_t          n_layers,
                  const fvm_coord_t         extrusion_vectors[],
                  const fvm_coord_t         distribution[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_NODAL_EXTRUDE_H__ */
