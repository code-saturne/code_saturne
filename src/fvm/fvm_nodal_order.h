#ifndef __FVM_NODAL_ORDER_H__
#define __FVM_NODAL_ORDER_H__

/*============================================================================
 * Ordering of nodal mesh entity lists and connectivity
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
 * Locally order cells and associated connectivity for a nodal mesh
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure.
 *   parent_global_number <-- global numbers of parent cells (if NULL, a
 *                            default 1 to n numbering is considered).
 *----------------------------------------------------------------------------*/

void
fvm_nodal_order_cells(fvm_nodal_t       *this_nodal,
                      const cs_gnum_t    parent_global_number[]);

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
                      const cs_gnum_t    parent_global_number[]);

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
                         const cs_gnum_t    parent_global_number[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_ORDER_H__ */
