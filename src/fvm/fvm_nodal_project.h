#ifndef __FVM_NODAL_PROJECT_H__
#define __FVM_NODAL_PROJECT_H__

/*============================================================================
 * Projection of nodal sections associated with a mesh
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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Project an extruded mesh to its base plane.
 *
 * parameters:
 *   this_nodal  <-> pointer to structure that should be cut in edges
 *   chosen_axis <-- indicate which axis is selected to extract edges
 *   error_count --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_project(fvm_nodal_t  *this_nodal,
                  int           chosen_axis,
                  cs_lnum_t    *error_count);

/*----------------------------------------------------------------------------
 * Reduce the spatial dimension of a mesh, discarding the last coordinate
 * component.
 *
 * The mesh's spatial dimension is reduced by 1.
 *
 * parameters:
 *   this_nodal <-> pointer to structure that projected
 *   matrix     <-- projection matrix
 *                  3D -> 2D: (a11, a12, a13, a21, a22, a23)
 *                  2D -> 1D: (a11, a12)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_project_coords(fvm_nodal_t  *this_nodal,
                         double       matrix[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_PROJECT_H__ */
