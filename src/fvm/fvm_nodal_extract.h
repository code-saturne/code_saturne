#ifndef __FVM_NODAL_EXTRACT_H__
#define __FVM_NODAL_EXTRACT_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/


/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Copy global vertex numbers to an array.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   g_vtx_num  --> global vertex numbers (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_global_vertex_num(const fvm_nodal_t  *this_nodal,
                                cs_gnum_t          *g_vtx_num);

/*----------------------------------------------------------------------------
 * Copy global element numbers of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the global element numbers are continued from one section to the next,
 * so to the user, all is as if the sections were concatenated.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   element_type  <-- type of elements to deal with
 *   g_elt_num     <-> pointer to global_element_num array (pre-allocated)
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_global_element_num(const fvm_nodal_t  *this_nodal,
                                 fvm_element_t       element_type,
                                 cs_gnum_t          *g_elt_num);

/*----------------------------------------------------------------------------
 * Copy vertex coordinates to an array.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   vertex_coords  --> vertices coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_vertex_coords(const fvm_nodal_t  *this_nodal,
                            cs_interlace_t      interlace,
                            cs_coord_t         *vertex_coords);

/*----------------------------------------------------------------------------
 * Copy element centers to an array.
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   entity_dim     <-- dimension of entities we want to count (0 to 3)
 *   cell_centers   --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_element_centers(const fvm_nodal_t  *this_nodal,
                              cs_interlace_t      interlace,
                              int                 entity_dim,
                              cs_coord_t         *cell_centers);

/*----------------------------------------------------------------------------
 * Copy element -> vertex connectivity of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the connectivity spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * Return local connectivity for sections of a given element_type.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   element_type  <-- type of elements of the section to deal with
 *   connectivity  <-> pointer to connectvity (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_strided_connect(const fvm_nodal_t  *this_nodal,
                              fvm_element_t       element_type,
                              cs_lnum_t          *connectivity);

/*----------------------------------------------------------------------------
 * Build inverse vertex -> element connectivity.
 *
 * The user is responsible for freeing the returned arrays.
 * The size of element_index[] should be n_vertices + 1, where n_vertices
 * is the value returned by fvm_nodal_get_n_entities(this_nodal, entity_dim).
 * The size of element_id[] should be element_index[n_vertices].
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * Note also that both the vertex -> element index and connectivity arrays
 * returned use 0 to n numbering.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   entity_dim    <-- dimension of entities we want to count (1 to 3)
 *   element_index --> vertex -> element connectivity index (O to n-1)
 *   element_id    --> vertex -> element connectivity (0 to n-1)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_vertex_elements(const fvm_nodal_t   *this_nodal,
                              int                  entity_dim,
                              cs_lnum_t          **element_index,
                              cs_lnum_t          **element_id);

/*----------------------------------------------------------------------------
 * Compute extents of a nodal mesh representation
 *
 * parameters:
 *   this_nodal   <-- pointer to mesh representation structure
 *   tolerance    <-- addition to local extents of each element:
 *                    extent = base_extent * (1 + tolerance)
 *   extents      <-> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_extents(const fvm_nodal_t  *this_nodal,
                  double              tolerance,
                  double              extents[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_EXTRACT_H__ */
