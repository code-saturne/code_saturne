#ifndef __FVM_NODAL_ORDER_H__
#define __FVM_NODAL_ORDER_H__

/*============================================================================
 * Ordering of nodal mesh entity lists and connectivity
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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
 * Locally order cells and associated connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent cells (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvm_nodal_order_cells(fvm_nodal_t       *this_nodal,
                      const fvm_gnum_t   parent_global_number[]);

/*----------------------------------------------------------------------------
 * Locally order faces and associated connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent faces (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvm_nodal_order_faces(fvm_nodal_t       *this_nodal,
                      const fvm_gnum_t   parent_global_number[]);

/*----------------------------------------------------------------------------
 * Locally order vertices and update connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent vertices (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvm_nodal_order_vertices(fvm_nodal_t       *this_nodal,
                         const fvm_gnum_t   parent_global_number[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_NODAL_ORDER_H__ */
