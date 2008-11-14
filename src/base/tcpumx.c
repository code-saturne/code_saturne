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

#undef _POSIX_SOURCE /* Otherwise compilation problem on VPP 5000 */
#undef _XOPEN_SOURCE /* Otherwise, compilation problem on SunOS */

/*============================================================================
 * Query time allocated to this process (useful mainly under PBS)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "tcpumx.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query CPU time allocated to this process
 *
 * Fortran interface:
 *
 * SUBROUTINE TCPUMX (TPS   , RET)
 * *****************
 *
 * DOUBLE PRECISION TPS        : <-- : remaining time (default: 7 days)
 * INTEGER          RET        : <-- : return code:
 *                             :     :  -1: error
 *                             :     :   0: no limit using this method
 *                             :     :   1: CPU limit determined
 *----------------------------------------------------------------------------*/

void CS_PROCF (tcpumx, TCPUMX) (double  *tps,
                                int     *ret)

{
  char * cs_maxtime;
  int    hrs, min, sec;
  int    nchamps = 0;

  *tps = 3600.0 * 24.0 * 7; /* "unlimited" values by default */
  *ret = 0;

  /* Get environment variable; for example, 100:10:10 */

  if ((cs_maxtime = getenv("CS_MAXTIME")) != NULL) {;

    nchamps = sscanf (cs_maxtime,"%d:%d:%d",&hrs,&min,&sec);

    /* If we only have 2 fields, they are hours and minutes (under PBS);
     * otherwise, if we do not have 3 fields, the information is unusable */

    if (nchamps == 2) {
      sec = 0;
      nchamps = 3;
    }

    /* Compute allocated CPU time in seconds */
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
