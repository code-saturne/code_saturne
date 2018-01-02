/*============================================================================
 * Resource allocation management (available time).
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#undef _POSIX_SOURCE /* Otherwise compilation problem on VPP 5000 */
#undef _XOPEN_SOURCE /* Otherwise, compilation problem on SunOS */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_SYS_RESOURCE_H)
#include <sys/resource.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_resource.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute remaining time allocated to this process when available
 * (usually, architectures other than IBM Blue Gene or Cray XT)
 *
 * parameters:
 *   tps <-> remaining time (default: 7 days)
 *
 * returns:
 *   -1 (error), 0 (no limit using this method), or 1 (limit determined)
 *----------------------------------------------------------------------------*/

static int
_t_remain(double  *tps)
{
  int retval = 0;

#if defined(HAVE_SYS_RESOURCE_H)

  struct rlimit ressources;
  struct rusage buf_time;
  struct rusage buf_time1;

  *tps = 3600.0 * 24.0 * 7; /* "unlimited" values by default */

  if ((retval = getrusage(RUSAGE_SELF, &buf_time)) < 0)
    bft_error(__FILE__, __LINE__, errno,
              "getrusage(RUSAGE_SELF) error.");

#if !defined(__bg__)
  else if ((retval = getrusage(RUSAGE_CHILDREN, &buf_time1)) < 0)
    bft_error(__FILE__, __LINE__, errno,
              "getrusage(RUSAGE_CHILDREN) error.");
#endif

  else if ((retval = getrlimit(RLIMIT_CPU, &ressources)) < 0)
    bft_error(__FILE__, __LINE__, errno,
              "getrlimit(RLIMIT_CPU) error.");

  /* If no error encountered (most probable case), use CPU limit returned by
     getrlimit (works at least with LSF batch systems under OSF1 or Linux),
     compute true remaining time, and put return code to 1 to indicate
     thet the remaining time is indeed limited. */

  if (retval == 0 && ressources.rlim_cur != RLIM_INFINITY) {
    *tps = (double)((int)ressources.rlim_cur
                    - (  buf_time.ru_utime.tv_sec  + buf_time.ru_stime.tv_sec
                       + buf_time1.ru_utime.tv_sec + buf_time1.ru_stime.tv_sec));
    retval = 1;
  }

#else

  retval = -1; /* getrusage(RUSAGE_SELF, ...) and getrlimit(RLIMIT_CPU, ...)
                  not available on this architecture */

#endif

  return retval;
}

/*----------------------------------------------------------------------------
 * Query CPU time allocated to this process
 *
 * parameters:
 *   tps <-> remaining time (default: 7 days)
 *
 * returns:
 *   -1 (error), 0 (no limit using this method), or 1 (limit determined)
 *----------------------------------------------------------------------------*/

static int
_t_cpu_max(double  *tps)
{
  char * cs_maxtime;
  int    hrs, min, sec;
  int    n_fields = 0;

  int retval = 0;

  *tps = 3600.0 * 24.0 * 7; /* "unlimited" values by default */

  /* Get environment variable; for example, 100:10:10 */

  if ((cs_maxtime = getenv("CS_MAXTIME")) != NULL) {;

    n_fields = sscanf (cs_maxtime,"%d:%d:%d",&hrs,&min,&sec);

    /* If we only have 2 fields, they are hours and minutes (under PBS);
     * otherwise, if we do not have 3 fields, the information is unusable */

    if (n_fields == 2) {
      sec = 0;
      n_fields = 3;
    }

    /* Compute allocated CPU time in seconds */
    if (n_fields == 3) {
      *tps = ((double)hrs)*3600. + ((double)min)*60. + ((double)sec);
      retval = 1;
#if 0
      printf("tcpumx n_fields = %d, hrs = %d, min = %d, sec = %d\n tps = %f\n",
             ret, hrs, min, sec, *tps);
#endif
    }
    else
      retval = -1;

  }

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Limit number of remaining time steps if the remaining allocated time is
 * too small to attain the requested number of steps.
 *
 * Fortran interface:
 *
 * subroutine armtsp (ntcabs, ntmabs)
 * *****************
 *
 * integer          ntcabs      : <-- : current time step number
 * integer          ntmabs      : <-> : maximum time step number
 *----------------------------------------------------------------------------*/

void CS_PROCF (armtps, ARMTPS)
(
 const cs_int_t  *ntcabs,
       cs_int_t  *ntmabs
)
{
  int ts_max = *ntmabs;

  cs_resource_get_max_timestep(*ntcabs, &ts_max);

  *ntmabs = ts_max;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Limit number of remaining time steps if the remaining allocated time is
 * too small to attain the requested number of steps.
 *
 * parameters:
 *   ts_cur <-- current time step number
 *   ts_max <-> maximum time step number
 *----------------------------------------------------------------------------*/

void
cs_resource_get_max_timestep(int   ts_cur,
                             int  *ts_max)
{
  /* Local variables */

  int t_lim_flag;
  double trestc, t_it_mean, alpha;
  double tcpuco, tresmn, titsmx;
  double t_margin = -1., tmoy00 = -1., t_it_prev = -1., t_it_sup = -1.;

  static int r_time_method = -1, ntcab0 = -1;
  static double trest0 = -1., tcpupr = -1.;

  if (ts_cur == *ts_max)
    return;

  /* Initialization at first pass */

  if (r_time_method ==  -1) {

    /* In parallel, rank 0 decides, broadcasts to others */

    if (cs_glob_rank_id <= 0) {

      /* First, try _t_remain. If no resource limits are set, we then try
         _t_cpu_max, which is based on the presence of the CS_MAXTIME
         environment variable */

      t_lim_flag = _t_remain(&trest0);
      if (t_lim_flag == 1)
        r_time_method = 1;
      else {
        t_lim_flag = _t_cpu_max(&trest0);
        if (t_lim_flag == 1)
          r_time_method = 2;
      }

    }

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      MPI_Bcast(&r_time_method, 1, MPI_INT, 0, cs_glob_mpi_comm);
#endif

    /* Remaining time and Wall-clock time at current iteration
       (which will become the previous one) */

    if (r_time_method > 0) {
      ntcab0 = ts_cur;
      tcpupr = cs_timer_wtime();
    }

  }

  /* At subsequent passes, use resource management method previously
     determined. */

  else if (r_time_method > 0) {

    /* In parallel, rank 0 decides, broadcasts to others */

    if (cs_glob_rank_id <= 0) {

      /* Mean time per iteration */

      /* previous iteration */
      tcpuco = cs_timer_wtime();
      t_it_prev = tcpuco - tcpupr;

      /* Current remaining time and mean iteration time */

      if (r_time_method == 1) {
        t_lim_flag = _t_remain(&trestc);
        tmoy00 = (trest0-trestc)/((double)(ts_cur-ntcab0));
      }
      else { /* if (r_time_method == 2) */
        /* Use initially allocated time */
        trestc = CS_MAX((trest0 - tcpuco), 0.);
        tmoy00 = tcpuco/((double)(ts_cur-ntcab0));
      }

      /* Estimate time per iteration (alpha > 0 safer) */

      alpha = 0.25;
      t_it_mean = alpha*t_it_prev + (1.-alpha)*tmoy00;

      /* Remaining time and CPU time at current iteration
         (which will become the previous one) */

      tcpupr = tcpuco;

      /* Time required for an additional iteration.
       *
       * Margin for I/O:
       * 100 times an iteration or 10% of time allocated to process
       * and at least 50 seconds or 1% of time allocated to process. */

      t_margin = CS_MIN(t_it_mean*100, trest0*0.1);
      t_margin = CS_MAX(t_margin, 50.);
      t_margin = CS_MAX(t_margin, trest0*0.01);

      /* Time for an additional iteration */

      t_it_sup = t_it_mean + t_margin;

      /* Test for calculation stop (in parallel, rank 0 decides). */

      tresmn = trestc;
      titsmx = t_it_sup;

      if (tresmn < titsmx) {
        *ts_max = ts_cur;
        bft_printf
          (_("===========================================================\n"
             "   ** Stop to avoid exceeding time allocation.\n"
             "      ----------------------------------------\n"
             "      maximum time step number set to: %d\n"
             "===========================================================\n"),
           *ts_max);
      }

    }

    /* Broadcast */

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      MPI_Bcast(ts_max, 1, CS_MPI_INT, 0, cs_glob_mpi_comm);
#endif

    if (cs_glob_rank_id <= 0 && ts_cur == *ts_max)
      bft_printf
        (_("===============================================================\n"
           "   ** Remaining time management\n"
           "      -------------------------\n"
           "      Remaining time allocated to the job       : ', %14.5e\n"
           "      Estimated time for another time step      : ', %14.5e\n"
           "        mean time for a time step               : ', %14.5e\n"
           "        time for the previous time step         : ', %14.5e\n"
           "        security margin                         : ', %14.5e\n"
           "===============================================================\n"),
         trestc, t_it_sup, tmoy00, t_it_prev, t_margin);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
