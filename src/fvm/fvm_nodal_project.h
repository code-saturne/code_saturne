#ifndef __FVM_NODAL_PROJECT_H__
#define __FVM_NODAL_PROJECT_H__

/*============================================================================
 * Projection of nodal sections associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2008  EDF

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
                  fvm_lnum_t   *error_count);

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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_NODAL_PROJECT_H__ */
