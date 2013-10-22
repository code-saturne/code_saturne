/*============================================================================
 * Utility functions for the diphasic lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_utils.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

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
 * Returns:
 *   1 if the two vertices are identical otherwise 0
 *----------------------------------------------------------------------------*/

int
cs_lagrang_check_colocalization(const cs_real_t p[3],
                                const cs_real_t q[3])
{
  int  sign = -1;

  /* We check if the two vertices are identical */

  if (   fabs(p[0] - q[0]) < 1e-15
      && fabs(p[1] - q[1]) < 1e-15
      && fabs(p[2] - q[2]) < 1e-15)
    sign = 1;
  else
    sign = 0;

  return  sign;
}

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
 *  an indicator on the orientation of the tetrahedron [p1, p2, p3, p4]
 *----------------------------------------------------------------------------*/

int
cs_lagrang_tetra_orientation(const cs_real_t  p1[3],
                             const cs_real_t  p2[3],
                             const cs_real_t  p3[3],
                             const cs_real_t  p4[3])
{
  int  returned_sign = -2; /* initialize to an incoherent value */

  /* points are assumed to be distinct */

  /*  | A B C |   | bx-ax  by-ay  bz-az |  */
  /*  | D E F | = | cx-ax  cy-ay  cz-az |  */
  /*  | G H I |   | dx-ax  dy-ay  dz-az |  */

  cs_real_t A = p2[0] - p1[0];
  cs_real_t B = p2[1] - p1[1];
  cs_real_t C = p2[2] - p1[2];
  cs_real_t D = p3[0] - p1[0];
  cs_real_t E = p3[1] - p1[1];
  cs_real_t F = p3[2] - p1[2];
  cs_real_t G = p4[0] - p1[0];
  cs_real_t H = p4[1] - p1[1];
  cs_real_t I = p4[2] - p1[2];

  /* minors computation */

  /*  M1 = D*H - E*G  */
  /*  M2 = A*H - B*G  */
  /*  M3 = A*E - B*D  */

  cs_real_t M1 = D*H - E*G;
  cs_real_t M2 = A*H - B*G;
  cs_real_t M3 = A*E - B*D;

  /* determinant computation */

  /*  det = (C*M1 - F*M2) + I*M3  */

  cs_real_t det = C*M1 - F*M2 + I*M3;

  /* *sign ==  0 => COPLANAR,
     *sign ==  1 => POSITIVE orientation,
     *sign == -1 => NEGATIVE orientation */

  if (det > 0.e0)
  {
    returned_sign = 1;
  } else if (det < 0.e0) {
    returned_sign = -1;
  } else {
    returned_sign = 0;
  }

  return  returned_sign;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
