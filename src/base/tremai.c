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

/*============================================================================
 * Compute remaining time allocated to this process
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

#undef _POSIX_SOURCE /* Otherwise compilation problem on VPP 5000 */
#undef _XOPEN_SOURCE /* Otherwise, compilation problem on SunOS */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "tremai.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute remaining time allocated to this process
 *
 * Fortran interface:
 *
 * SUBROUTINE TREMAI (TPS   , RET)
 * *****************
 *
 * DOUBLE PRECISION TPS        : <-- : remaining time (default: 7 days)
 * INTEGER          RET        : <-- : return code:
 *                             :     :  -1: error
 *                             :     :   0: no limit using this method
 *                             :     :   1: CPU limit determined
 *----------------------------------------------------------------------------*/

void CS_PROCF (tremai, TREMAI) (double  *tps,
                                int     *ret)
{
  struct rlimit ressources;
  struct rusage buf_time;
  struct rusage buf_time1;

  *tps = 3600.0 * 24.0 * 7; /* "unlimited" values by default */

/* Architectures other than IBM Blue Gene or Cray XT */
#if   !defined(__blrts__) && !defined(__bgp__) \
   && !defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  if ((*ret = getrusage(RUSAGE_SELF, &buf_time)) < 0)
    fprintf(stderr, "getrusage(RUSAGE_SELF) error:\n%s\n", strerror(errno));
  else if ((*ret = getrusage(RUSAGE_CHILDREN, &buf_time1)) < 0)
    fprintf(stderr, "getrusage(RUSAGE_CHILDREN) error:\n%s\n", strerror(errno));
  else if ((*ret = getrlimit(RLIMIT_CPU, &ressources)) < 0)
    fprintf(stderr, "getrlimit(RLIMIT_CPU) error:\n%s\n", strerror(errno));

  /* If no error encountered (most probable case), use CPU limit returned by
     getrlimit (works at least with LSF batch systems under OSF1 or Linux),
     compute true remaining time, and put return code to 1 to indicate
     thet the remaining time is indeed limited. */

  if (*ret == 0 && ressources.rlim_cur != RLIM_INFINITY) {
    *tps = (double)((int)ressources.rlim_cur
                    - (  buf_time.ru_utime.tv_sec  + buf_time.ru_stime.tv_sec
                       + buf_time1.ru_utime.tv_sec + buf_time1.ru_stime.tv_sec));
    *ret = 1;
  }

#else /* IBM Blue Gene or Cray XT */

  *ret = -1; /* getrusage(RUSAGE_SELF, ...) and getrlimit(RLIMIT_CPU, ...)
                not available on this architecture */

#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
