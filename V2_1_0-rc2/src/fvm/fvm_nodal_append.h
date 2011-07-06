#ifndef __FVM_NODAL_APPEND_H__
#define __FVM_NODAL_APPEND_H__

/*============================================================================
 * Append sections to a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2006  EDF

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
 * Append a new section to an existing fvm_nodal mesh, and transfer
 * ownership of the given connectivity and optional parent number arrays to
 * that section.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvm_nodal_append_by_transfer(fvm_nodal_t    *this_nodal,
                             fvm_lnum_t      n_elements,
                             fvm_element_t   type,
                             fvm_lnum_t      face_index[],
                             fvm_lnum_t      face_num[],
                             fvm_lnum_t      vertex_index[],
                             fvm_lnum_t      vertex_num[],
                             fvm_lnum_t      parent_element_num[]);

/*----------------------------------------------------------------------------
 * Append a new section to an existing fvm_nodal mesh, sharing the given
 * given connectivity and optional parent number arrays with the caller.
 *
 * The caller should not destroy or modify the arrays passed to this
 * function until the nodal mesh is destroyed.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvm_nodal_append_shared(fvm_nodal_t    *this_nodal,
                        fvm_lnum_t      n_elements,
                        fvm_element_t   type,
                        fvm_lnum_t      face_index[],
                        fvm_lnum_t      face_num[],
                        fvm_lnum_t      vertex_index[],
                        fvm_lnum_t      vertex_num[],
                        fvm_lnum_t      parent_element_num[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_NODAL_APPEND_H__ */
