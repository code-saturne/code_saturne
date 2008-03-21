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

/*#undef _POSIX_SOURCE / * Sinon, problème de compilation sur VPP 5000 * /
#undef _XOPEN_SOURCE / * Sinon, problème de compilation sur SunOS    */


/* includes système */

#include <stdio.h>
#include <stdlib.h>

/* Includes librairie */

#include "cs_base.h"
#include "tcpumx.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 *  Récuperation du temps cpu alloué au process
 *  (utile notamment pour cluster sous PBS)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Récuperation du temps cpu alloué au process
 *
 * Interface Fortran :
 *
 * SUBROUTINE TCPUMX (TPS   , RET)
 * *****************
 *
 * DOUBLE PRECISION TPS        : <-- : Temps restant (défaut : 7 jours)
 * INTEGER          RET        : <-- : Code de retour ;
 *                             :     :  -1 : erreur
 *                             :     :   0 : pas de limite via cette méthode
 *                             :     :   1 : limite de temps CPU déterminée
 *----------------------------------------------------------------------------*/

void CS_PROCF (tcpumx, TCPUMX) (double  *tps,
                                int     *ret)

{
  char * cs_maxtime;
  int    hrs, min, sec;
  int    nchamps = 0;

  *tps = 3600.0 * 24.0 * 7; /* valeur "illimitée" par défaut */
  *ret = 0;

  /* Récuperation de la variable d'environnement ; ex : 100:10:10 */

  if ((cs_maxtime = getenv("CS_MAXTIME")) != NULL) {;

    nchamps = sscanf (cs_maxtime,"%d:%d:%d",&hrs,&min,&sec);

    /* Si on n'a que 2 champs ce sont les heures et les minutes (avec PBS) ;
     * sinon si l'on n'a pas 3 champs l'information n'est pas exploitable */

    if (nchamps == 2) {
      sec = 0;
      nchamps = 3;
    }

    /* Calcul du temps CPU alloué en secondes */
    if (nchamps == 3) {
      *tps = ((double)hrs)*3600. + ((double)min)*60. + ((double)sec);
      *ret = 1;
#if 0
      printf("tcpumx nchamps = %d,hrs = %d, min = %d, sec = %d\n tps = %f\n",
             ret, hrs, min, sec, *tps);
#endif
    }
    else
      *ret = -1;

  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

