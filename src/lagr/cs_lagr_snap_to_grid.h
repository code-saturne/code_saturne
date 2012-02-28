#ifndef __CS_LAGR_SNTOGR_H__
#define __CS_LAGR_SNTOGR_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if IEEE 754 standard is respected for floating storage for this
 * architecture. If the standard is not respected the particle trajectography
 * may be wrong.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (csieee,CSIEEE)(void);

/*----------------------------------------------------------------------------
 * Check the relative localization of two vertices. We want to know if these
 * two vertices are identical.
 *
 * pvalmax              --> upperbound on coordinates
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
CS_PROCF (colosg,COLOSG)(cs_real_t  *pvalmax,
                         cs_real_t  *px,
                         cs_real_t  *py,
                         cs_real_t  *pz,
                         cs_real_t  *qx,
                         cs_real_t  *qy,
                         cs_real_t  *qz,
                         cs_int_t   *sign);

/*----------------------------------------------------------------------------
 * Look for coordinate system orientation to locate particles in relation to
 * faces.
 *
 * pvalmax              --> upper bound on coordinates
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
 * pturb                <->
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (cotusg,COTUSG)(cs_real_t   *pvalmax,
                         cs_real_t   *px,
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
                         cs_int_t    *sign,
                         cs_int_t    *pturb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_SNTOGR_H__ */
