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

#ifndef __CS_TCPUMX_H__
#define __CS_TCPUMX_H__

/* Includes librairie */

#include "cs_base.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 *  Recuperation du temps cpu alloue au process
 *  (utile notamment pour cluster sous PBS)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Recuperation du temps cpu alloue au process
 *
 * Interface Fortran :
 *
 * SUBROUTINE TCPUMX (TPS   , RET)
 * *****************
 *
 * DOUBLE PRECISION TPS        : <-- : Temps restant (defaut : 7 jours)
 * INTEGER          RET        : <-- : Code de retour ;
 *                             :     :  -1 : erreur
 *                             :     :   0 : pas de limite via cette methode
 *                             :     :   1 : limite de temps CPU determinee
 *----------------------------------------------------------------------------*/

void CS_PROCF (tcpumx, TCPUMX) (double  *tps,
                                int     *ret);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_TCPUMX_H_ */
