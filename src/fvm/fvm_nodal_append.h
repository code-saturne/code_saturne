#ifndef __FVM_NODAL_APPEND_H__
#define __FVM_NODAL_APPEND_H__

/*============================================================================
 * Append sections to a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
                             cs_lnum_t       n_elements,
                             fvm_element_t   type,
                             cs_lnum_t       face_index[],
                             cs_lnum_t       face_num[],
                             cs_lnum_t       vertex_index[],
                             cs_lnum_t       vertex_num[],
                             cs_lnum_t       parent_element_num[]);

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
                        cs_lnum_t       n_elements,
                        fvm_element_t   type,
                        cs_lnum_t       face_index[],
                        cs_lnum_t       face_num[],
                        cs_lnum_t       vertex_index[],
                        cs_lnum_t       vertex_num[],
                        cs_lnum_t       parent_element_num[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_APPEND_H__ */
