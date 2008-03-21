/*============================================================================
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

#ifndef __CS_TREMAI_H__
#define __CS_TREMAI_H__


/* Includes librairie */

#include "cs_base.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 *  Prototype de la fonction tremai permettant de connaitre le temps restant
 *  alloue au process
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Calcul du temps restant alloué au process
 *
 * Interface Fortran :
 *
 * SUBROUTINE TREMAI (TPS   , RET)
 * *****************
 *
 * DOUBLE PRECISION TPS        : <-- : Temps restant (défaut : 7 jours)
 * INTEGER          RET        : <-- : Code de retour ;
 *                             :     :  -1 : erreur
 *                             :     :   0 : pas de limite via cette méthode
 *                             :     :   1 : limite de temps CPU déterminée
 *----------------------------------------------------------------------------*/

void CS_PROCF (tremai, TREMAI) (double  *tps,
                                int     *ret);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_TREMAI_H_ */
