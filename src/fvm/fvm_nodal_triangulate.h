#ifndef __FVM_NODAL_TRIANGULATE_H__
#define __FVM_NODAL_TRIANGULATE_H__

/*============================================================================
 * Triangulation of nodal sections associated with a mesh
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * Triangulate all sections of a nodal mesh.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be triangulated
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_triangulate(fvm_nodal_t  *this_nodal,
                      cs_lnum_t    *error_count);

/*----------------------------------------------------------------------------
 * Triangulate polygonal sections of a nodal mesh.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be triangulated
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_triangulate_polygons(fvm_nodal_t  *this_nodal,
                               cs_lnum_t    *error_count);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_NODAL_TRIANGULATE_H__ */
