/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_LAGRANG_H__
#define __CS_LAGRANG_H__

/*============================================================================
 * Utilitarian functions for the diphasic lagrangian module
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
CS_PROCF (coloca,COLOCA)(cs_real_t  *pvalmax,
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
CS_PROCF (coturn,COTURN)(cs_real_t   *pvalmax,
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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_LAGRANG_H__ */
