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

#ifndef __CS_SELECTOR_H__
#define __CS_SELECTOR_H__

/*============================================================================
 * Structure principale associée à un maillage
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfbr, CSGFBR)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const faces_number,   /* --> faces number */
 int                 *const faces         /* --> faces  */
 CS_ARGF_SUPP_CHAINE
 );

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfac, CSGFAC)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const faces_number,   /* --> faces number */
 int                 *const faces         /* --> faces  */
 CS_ARGF_SUPP_CHAINE
 );

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgcel, CSGCEL)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const cells_number,   /* --> cells number */
 int                 *const cells         /* --> cells  */
 CS_ARGF_SUPP_CHAINE
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SELECTOR_H__ */
