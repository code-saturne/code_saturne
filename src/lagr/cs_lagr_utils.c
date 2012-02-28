/*============================================================================
 * Utilitarian functions for the diphasic lagrangian module
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

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
 * px                   --> X coordinate of the vertex P
 * py                   --> Y coordinate of the vertex P
 * pz                   --> Z coordinate of the vertex P
 * qx                   --> X coordinate of the vertex Q
 * qy                   --> Y coordinate of the vertex Q
 * qz                   --> Z coordinate of the vertex Q
 * sign                 <-> return tag (1 -> identical else 0)
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (coloca,COLOCA)(cs_real_t  *px,
                         cs_real_t  *py,
                         cs_real_t  *pz,
                         cs_real_t  *qx,
                         cs_real_t  *qy,
                         cs_real_t  *qz,
                         cs_int_t   *sign)
{


  /* We check if the two vertices are identical */



  if (fabs(px-qx) < 1e-15 && fabs(py-qy) < 1e-15 && fabs(pz-qz) < 1e-15)
    *sign = 1;
  else
    *sign = 0;
}

/*----------------------------------------------------------------------------
 * Look for coordinate system orientation to locate particles in relation to
 * faces.
 *
 * px                   --> X coordinate of the first vertex
 * py                   --> Y coordinate of the first vertex
 * pz                   --> Z coordinate of the first vertex
 * qx                   --> X coordinate of the second vertex
 * qy                   --> Y coordinate of the second vertex
 * qz                   --> Z coordinate of the second vertex
 * cdgx                 --> X coordinate of the third vertex
 * cdgy                 --> Y coordinate of the third vertex
 * cdgz                 --> Z coordinate of the third vertex
 * crgx                 --> X coordinate of the fourth vertex
 * crgy                 --> Y coordinate of the fourth vertex
 * crgz                 --> Z coordinate of the fourth vertex
 * sign                 <-> orientation of the four vertices.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (coturn,COTURN)(cs_real_t   *px,
                         cs_real_t   *py,
                         cs_real_t   *pz,
                         cs_real_t   *qx,
                         cs_real_t   *qy,
                         cs_real_t   *qz,
                         cs_real_t   *cdgx,
                         cs_real_t   *cdgy,
                         cs_real_t   *cdgz,
                         cs_real_t   *crdx,
                         cs_real_t   *crdy,
                         cs_real_t   *crdz,
                         cs_int_t    *sign)
{

  cs_real_t ax = *px;
  cs_real_t ay = *py;
  cs_real_t az = *pz;
  cs_real_t bx = *qx;
  cs_real_t by = *qy;
  cs_real_t bz = *qz;
  cs_real_t cx = *cdgx;
  cs_real_t cy = *cdgy;
  cs_real_t cz = *cdgz;
  cs_real_t dx = *crdx;
  cs_real_t dy = *crdy;
  cs_real_t dz = *crdz;


  /* points are assumed to be distincts */

  /*  | A B C |   | bx-ax  by-ay  bz-az |  */
  /*  | D E F | = | cx-ax  cy-ay  cz-az |  */
  /*  | G H I |   | dx-ax  dy-ay  dz-az |  */

  cs_real_t A = bx - ax;
  cs_real_t B = by - ay;
  cs_real_t C = bz - az;
  cs_real_t D = cx - ax;
  cs_real_t E = cy - ay;
  cs_real_t F = cz - az;
  cs_real_t G = dx - ax;
  cs_real_t H = dy - ay;
  cs_real_t I = dz - az;


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
    *sign = 1;
  } else if (det < 0.e0) {
    *sign = -1;
  } else {
    *sign = 0;
  }


}

/*----------------------------------------------------------------------------*/

END_C_DECLS
