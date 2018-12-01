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

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_resource.c
 *
 *  \brief Resource allocation management (available time).
 */

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

static double  _wt_limit = -1.;
static double  _wt_safe = 0.95;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute remaining time allocated to this process when available.
 *
 * parameters:
 *   tps <-> remaining time (default: 30 days)
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

  *tps = 3600.0 * 24.0 * 30; /* "unlimited" values by default */

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
     that the remaining time is indeed limited. */

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
 * Initialize CPU time allocated based on CS_MAXTIME and user settings.
 *----------------------------------------------------------------------------*/

static void
_init_wt_limit(void)
{
  /* Get environment variable; for example, 100:10:10 */

  const char *cs_maxtime = getenv("CS_MAXTIME");

  if (cs_maxtime != NULL) {

    int hrs = -1, min = -1, sec = -1;
    int n_fields = sscanf(cs_maxtime,"%d:%d:%d",&hrs,&min,&sec);

    /* If we only have 1 field, assume it is in seconds (user defined) */

    if (n_fields == 1) {
      int n_secs = hrs;
      hrs = n_secs / 3600;
      n_secs = n_secs % 3600;
      min = n_secs / 60 ;
      sec = n_secs %60;
      n_fields = 3;
    }

    /* If we only have 2 fields, they are hours and minutes (under PBS);
     * otherwise, if we do not have 3 fields, the information is unusable */

    if (n_fields == 2) {
      sec = 0;
      n_fields = 3;
    }

    /* Compute allocated CPU time in seconds */
    if (n_fields == 3) {
      _wt_limit = ((double)hrs)*3600. + ((double)min)*60. + ((double)sec);
      bft_printf(_("\n Wall-clock time limit set by CS_MAXTIME: %dh:%dm:%ds\n"),
                 hrs, min, sec);
    }
    else {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("\n%s: Failed to parse CS_MAXTIME = \"%s\"\n"),
                 __func__, cs_maxtime);
    }

  }
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current wall-clock time limit.
 *
 * \return current wall-time limit (in seconds), or -1
 */
/*----------------------------------------------------------------------------*/

double
cs_resource_get_wt_limit(void)
{
  return _wt_limit;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set wall-clock time limit.
 *
 * \param[in]  wt  wall-time limit (in seconds), or -1
 */
/*----------------------------------------------------------------------------*/

void
cs_resource_set_wt_limit(double  wt)
{
  _wt_limit = wt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Limit number of remaining time steps if the remaining allocated
 *        time is too small to attain the requested number of steps.
 *
 * \param[in]       ts_cur  current time step number
 * \param[in, out]  ts_max  maximum time step number
 */
/*----------------------------------------------------------------------------*/

void
cs_resource_get_max_timestep(int   ts_cur,
                             int  *ts_max)
{
  /* Local variables */

  int t_lim_flag;

  static int r_time_method = -1;
  static double wtrem0 = -1.;

  if (ts_cur == *ts_max)
    return;

  /* Initialization at first pass */

  if (r_time_method ==  -1) {

    r_time_method = 0;

    /* In parallel, rank 0 decides, broadcasts to others */

    if (cs_glob_rank_id <= 0) {

      /* First, try _t_remain. If no resource limits are set, we then use
         _t_wt_limit, which is based on the presence of the CS_MAXTIME
         environment variable and on user settings. */

      t_lim_flag = _t_remain(&wtrem0);
      if (t_lim_flag == 1)
        r_time_method = 1;

      _init_wt_limit();

    }

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      MPI_Bcast(&r_time_method, 1, MPI_INT, 0, cs_glob_mpi_comm);
      MPI_Bcast(&_wt_limit, 1, MPI_DOUBLE, 0, cs_glob_mpi_comm);
    }
#endif

  }

  /* Use resource management method previously determined. */

  if (r_time_method > 0 || _wt_limit > 0) {

    /* In parallel, rank 0 decides, broadcasts to others */

    if (cs_glob_rank_id <= 0) {

      /* Current wall-clock and remaining time */

      double wt_cur = cs_timer_wtime();
      double wt_rem = -1;

      if (r_time_method == 1)
        t_lim_flag = _t_remain(&wt_rem);
      else if (r_time_method == 2) {
        /* Use initially allocated time */
        wt_rem = CS_MAX((wtrem0 - wt_cur), 0.);
      }

      if (_wt_limit > 0) {
        double wt_rem_l = _wt_limit - wt_cur;
        if (wt_rem_l < wt_rem || wt_rem < 0)
          wt_rem = wt_rem_l;
      }

      if (wt_cur >= (wt_rem + wt_cur)*_wt_safe) {

        *ts_max = ts_cur;

        bft_printf
          (_("===========================================================\n"
             "   ** Stop to avoid exceeding time allocation.\n"
             "      ----------------------------------------\n"
             "      maximum time step number set to: %d\n"
             "===========================================================\n"),
           *ts_max);

        FILE *_f = fopen("run_status.exceeded_time_limit", "w");
        if (_f != NULL) {
          fprintf(_f, "%d\n", ts_cur);
          fclose(_f);
        }

      }

    }

    /* Broadcast */

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      MPI_Bcast(ts_max, 1, CS_MPI_INT, 0, cs_glob_mpi_comm);
#endif

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
