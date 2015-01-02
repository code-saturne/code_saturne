#ifndef __CS_LAGR_UTILS_H__
#define __CS_LAGR_UTILS_H__

/*============================================================================
 * Utility functions for the diphasic lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check the relative localization of two vertices. We want to know if these
 * two vertices are identical.
 *
 * parameters:
 *   p             <-- X, Y, Z coordinates of the vertex P
 *   q             <-- X, Y, Z coordinates of the vertex Q
 *
 * returns:
 *   1 if the two vertices are identical otherwise 0
 *----------------------------------------------------------------------------*/

int
cs_lagrang_check_colocalization(const double  p[3],
                                const double  q[3]);

/*----------------------------------------------------------------------------
 * Look for coordinate system orientation to locate particles in relation to
 * faces.
 *
 * parameters:
 *   p1            <-- X, Y, Z coordinate of the first vertex
 *   p2            <-- X, Y, Z coordinate of the second vertex
 *   p3            <-- X, Y, Z coordinate of the third vertex
 *   p4            <-- X, Y, Z coordinate of the fourth vertex
 *
 * returns:
 *   an indicator on the orientation of the tetrahedron [p1, p2, p3, p4]
 *----------------------------------------------------------------------------*/

int
cs_lagrang_tetra_orientation(const double  p1[3],
                             const double  p2[3],
                             const double  p3[3],
                             const double  p4[3]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_UTILS_H__ */
