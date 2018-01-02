/*============================================================================
 * Program timing information
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

#include "ecs_def.h"

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */

#include <time.h>

#if defined (HAVE_GETTIMEOFDAY)
#include <sys/time.h>
#endif

#if defined (HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#elif defined(_POSIX_SOURCE)
#include <sys/times.h>
#include <unistd.h>
#endif

/*
 * Optional library and ECS headers
 */

#include "ecs_def.h"
#include "ecs_timer.h"

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static _Bool _ecs_timer_initialized = false;

/* Wall-clock time */

#if defined (HAVE_GETTIMEOFDAY)
static struct timeval  _ecs_timer_wtime_tv_start;
#else
static time_t _ecs_timer_wtime_start;
#endif

/* CPU time */

#if defined (HAVE_GETRUSAGE)
#elif defined(_POSIX_SOURCE)
static time_t _ecs_timer_unit = 0;
#else
static clock_t _ecs_timer_clock_start;
#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

static void
_ecs_timer_initialize(void)
{
#if defined (HAVE_GETTIMEOFDAY)
  (void)gettimeofday(&_ecs_timer_wtime_tv_start, NULL);
#else
  (void)time(&_ecs_timer_wtime_start);
#endif

#if defined (HAVE_GETRUSAGE)
#elif defined(_POSIX_SOURCE)
  _ecs_timer_unit = (double)sysconf(_SC_CLK_TCK);
#else
  _ecs_timer_clock_start = clock();
#endif

  _ecs_timer_initialized = true;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Return Wall clock time
 *
 * \return elapsed time from first call of a function of the ecs_timer_...()
 *         series, or -1 if unable to compute.
 */

double
ecs_timer_wtime(void)
{
  double this_wtime = -1.;

  /* Ensure initialization */

  if (_ecs_timer_initialized == false)
    _ecs_timer_initialize();

  /* Compute elapsed time */

#if defined (HAVE_GETTIMEOFDAY)

 {
    struct timeval  wtime_tv_current;

    if (gettimeofday(&wtime_tv_current, NULL) == 0) {

      /* Perform carry for later subtraction */
      if (_ecs_timer_wtime_tv_start.tv_usec > wtime_tv_current.tv_usec) {
        int nsec = (_ecs_timer_wtime_tv_start.tv_usec - wtime_tv_current.tv_usec)
                   / 1000000 + 1;
        wtime_tv_current.tv_usec += 1000000 * nsec;
        wtime_tv_current.tv_sec -= nsec;
      }
      if (  wtime_tv_current.tv_usec - _ecs_timer_wtime_tv_start.tv_usec
          > 1000000) {
        int nsec = (wtime_tv_current.tv_usec - _ecs_timer_wtime_tv_start.tv_usec)
                   / 1000000;
        wtime_tv_current.tv_usec -= 1000000 * nsec;
        wtime_tv_current.tv_sec += nsec;
      }

      this_wtime =   (  wtime_tv_current.tv_sec
                      - _ecs_timer_wtime_tv_start.tv_sec)
                   + (  wtime_tv_current.tv_usec
                      - _ecs_timer_wtime_tv_start.tv_usec) * 1.e-6 ;

    }

 }

#else

 {
   time_t wtime_current;

   if (time(&wtime_current) != (time_t)-1)
     this_wtime = difftime(wtime_current, _ecs_timer_wtime_start);
 }

#endif

  return this_wtime;
}

/*!
 * \brief Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see ecs_timer_cpu_time_method()), at least one of
 * the ecs_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * \return current CPU time usage, or -1 if unable to compute.
 */

double
ecs_timer_cpu_time(void)
{
  double cpu_time = -1.;

  /* Ensure initialization */

  if (_ecs_timer_initialized == 0)
    _ecs_timer_initialize();

  /* Compute CPU time */

#if defined (HAVE_GETRUSAGE)

 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     cpu_time  =    usage.ru_utime.tv_sec  + usage.ru_stime.tv_sec
                 + (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) * 1.e-6;
   }
 }

#elif defined(_POSIX_SOURCE)

 {
    static struct tms  ptimer;

    if (_ecs_timer_unit != -1 && times(&ptimer) != -1) {
      cpu_time =   ((double)(ptimer.tms_utime + ptimer.tms_stime))
                 / _ecs_timer_unit;
    }

 }

#else /* Use minimal C library function */

  if (_ecs_timer_clock_start != -1) {

    static clock_t  clock_current;

    clock_current = clock();
    if (clock_current != (clock_t)-1)
      cpu_time
        = ((double)(clock_current - _ecs_timer_clock_start)) / CLOCKS_PER_SEC;

  }

#endif

  return cpu_time;
}

/*!
 * \brief Return separate user and system CPU times.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available, this function will return -1 values.
 *
 * \param [out] user_time current user CPU usage.
 * \param [out] system_time current system CPU usage.
 */

void
ecs_timer_cpu_times(double *user_time,
                    double *system_time)
{
  /* Ensure initialization */

  if (_ecs_timer_initialized == 0)
     _ecs_timer_initialize();

  *user_time   = -1.;
  *system_time = -1.;

  /* Compute CPU time */

#if defined (HAVE_GETRUSAGE)

 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     *user_time    = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-6;
     *system_time  = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.e-6;
   }
 }

#elif defined(_POSIX_SOURCE)

 {
    static struct tms  ptimer;

    if (_ecs_timer_unit != -1 && times(&ptimer) != -1) {
      *user_time    = ((double)ptimer.tms_utime)  / _ecs_timer_unit;
      *system_time  = ((double)ptimer.tms_stime)  / _ecs_timer_unit;
    }

 }

#endif
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
