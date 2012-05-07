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

#undef _POSIX_SOURCE /* Sinon, problème de compilation sur VPP 5000 */
#undef _XOPEN_SOURCE /* Sinon, problème de compilation sur SunOS    */


/* includes système */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>


/* Includes librairie */

#include "cs_base.h"
#include "tremai.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 *  Calcul du temps restant alloué au process
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
                                int     *ret)
{
  struct rlimit ressources;
  struct rusage buf_time;
  struct rusage buf_time1;

  *tps = 3600.0 * 24.0 * 7; /* valeur "illimitée" par défaut */

/* Architectures hors IBM Blue Gene ou Cray XT */
#if   !defined(__blrts__) && !defined(__bg__) \
   && !defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  if ((*ret = getrusage(RUSAGE_SELF, &buf_time)) < 0)
    fprintf(stderr, "getrusage(RUSAGE_SELF) error:\n%s\n", strerror(errno));
  else if ((*ret = getrusage(RUSAGE_CHILDREN, &buf_time1)) < 0)
    fprintf(stderr, "getrusage(RUSAGE_CHILDREN) error:\n%s\n", strerror(errno));
  else if ((*ret = getrlimit(RLIMIT_CPU, &ressources)) < 0)
    fprintf(stderr, "getrlimit(RLIMIT_CPU) error:\n%s\n", strerror(errno));

  /* Si aucune erreur (le plus probable) et limitation de temps CPU
     indiquée par getrlimit (normalement le cas avec système de batch LSF
     sous OSF1 ou Linux par exemple), on calcule le temps restant réel,
     et on met le code retour à 1 pour indiquer que le temps restant
     est effectivement limité */

  if (*ret == 0 && ressources.rlim_cur != RLIM_INFINITY) {
    *tps = (double)((int)ressources.rlim_cur
                    - (  buf_time.ru_utime.tv_sec  + buf_time.ru_stime.tv_sec
                       + buf_time1.ru_utime.tv_sec + buf_time1.ru_stime.tv_sec));
    *ret = 1;
  }

#else /* IBM Blue Gene ou Cray XT */

  *ret = -1; /* getrusage(RUSAGE_SELF, ...) et getrlimit(RLIMIT_CPU, ...)
                non disponibles sur cette architecture */

#endif

#if defined(_CS_HAVE_MPI) /* Ensure all ranks have the same info
                             (especially for LoadLeveler) */
  if (cs_glob_base_nbr > 1) {
    double buf[2];
    buf[0] = *ret; buf[1] = *tps;
    MPI_Bcast(buf, 2, MPI_DOUBLE, 0, cs_glob_base_mpi_comm);
    *ret = buf[0]; *tps = buf[1];
  }
#endif

}

#ifdef __cplusplus
}
#endif /* __cplusplus */

