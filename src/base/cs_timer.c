/*============================================================================
 * Program timing information
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#if defined(HAVE_CONFIG_H)
#  include "cs_config.h"
#endif

#if defined(HAVE_CLOCK_GETTIME)
#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif
#endif

/*-----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <time.h>

#if defined (HAVE_GETTIMEOFDAY)
#include <sys/time.h>
#endif

#if defined (HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(_POSIX_SOURCE)
#include <sys/times.h>
#include <unistd.h>
#endif

/* Disable automatically-defined HAVE_CLOCK_GETTIME on Cygwin */

#if defined(HAVE_CLOCK_GETTIME) && defined(__CYGWIN__)
#undef HAVE_CLOCK_GETTIME
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Timing methods */

typedef enum {

  CS_TIMER_DISABLE,
  CS_TIMER_CLOCK_GETTIME,
  CS_TIMER_GETTIMEOFDAY,
  CS_TIMER_GETRUSAGE,
  CS_TIMER_TIME,
  CS_TIMER_TIMES,
  CS_TIMER_CLOCK

} cs_timer_method_t;

/* Function pointer types */

typedef void
(_cs_timer_wall_func_t) (cs_timer_t  *timer);

typedef void
(_cs_timer_cpu_func_t) (cs_timer_t  *timer);

/*-------------------------------------------------------------------------------
 * Local macro documentation
 *-----------------------------------------------------------------------------*/

/*! \fn CS_TIMER_COUNTER_INIT(_t)
 *
 * \brief Initialize timer counter.
 *
 * \param [out] _t  resulting counter.
 */

/*! \fn CS_TIMER_COUNTER_ADD(_res, _c0, _c1)
 *
 * \brief Add timer counter.
 *
 * The result may be identical to one of the 2 counters to add.
 *
 * \param [out] _res  resulting counter.
 * \param [in]  _c0   counter to add.
 * \param [in]  _c1   counter to add.
 */

/*-----------------------------------------------------------------------------
 * Local function prototypes (descriptions later).
 *-----------------------------------------------------------------------------*/

/* Wall-clock timer functions */

void
_cs_timer_wall_null(cs_timer_t  *timer);

#if defined (HAVE_CLOCK_GETTIME)
void
_cs_timer_wall_clock_gettime(cs_timer_t  *timer);
#endif

#if defined (HAVE_GETTIMEOFDAY)
void
_cs_timer_wall_gettimeofday(cs_timer_t  *timer);
#endif

void
_cs_timer_wall_stdc_time(cs_timer_t  *timer);

/* CPU timer functions */

void
_cs_timer_cpu_null(cs_timer_t  *timer);

#if defined (HAVE_CLOCK_GETTIME)
void
_cs_timer_cpu_clock_gettime(cs_timer_t  *timer);
#endif

#if defined (HAVE_GETRUSAGE)
void
_cs_timer_cpu_getrusage(cs_timer_t  *timer);
#endif

#if defined(_POSIX_SOURCE)
void
_cs_timer_cpu_times(cs_timer_t  *timer);
#endif

void
_cs_timer_cpu_stdc_clock(cs_timer_t  *timer);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static _Bool _cs_timer_initialized = false;
static cs_timer_method_t _cs_timer_wall_method = CS_TIMER_DISABLE;
static cs_timer_method_t _cs_timer_cpu_method = CS_TIMER_DISABLE;
static _cs_timer_wall_func_t *_cs_timer_wall = _cs_timer_wall_null;
static _cs_timer_cpu_func_t *_cs_timer_cpu = _cs_timer_cpu_null;

/* Wall-clock time */

static time_t _cs_timer_stdc_time_start;

/* CPU time */

#if defined(_POSIX_SOURCE)
static long _cs_timer_unit = 0;
#endif

static clock_t _cs_timer_clock_start;

/* Reference times */

static cs_timer_t  _cs_timer_start;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Null function to get wall-clock time.
 *
 * parameters:
 *   timer  <-- pointer to timing structure (ignored)
 *----------------------------------------------------------------------------*/

void
_cs_timer_wall_null(cs_timer_t  *timer)
{
  CS_UNUSED(timer);
}

#if defined (HAVE_CLOCK_GETTIME)

/*----------------------------------------------------------------------------
 * Get wall-clock time using clock_gettime().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_wall_clock_gettime(cs_timer_t  *timer)
{
  struct timespec w_time;
  (void)clock_gettime(CLOCK_REALTIME, &w_time);
  timer->wall_sec = w_time.tv_sec;
  timer->wall_nsec = w_time.tv_nsec;
}

#endif /* defined (HAVE_CLOCK_GETTIME) */

#if defined (HAVE_GETTIMEOFDAY)

/*----------------------------------------------------------------------------
 * Get wall-clock time using gettimeofday().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_wall_gettimeofday(cs_timer_t  *timer)
{
  struct timeval  tv_time;
  (void)gettimeofday(&tv_time, NULL);
  timer->wall_sec = tv_time.tv_sec;
  timer->wall_nsec = tv_time.tv_usec*1000;
}

#endif /* defined (HAVE_GETTIMEOFDAY) */

/*----------------------------------------------------------------------------
 * Get wall-clock time using time().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_wall_stdc_time(cs_timer_t  *timer)
{
  time_t wtime_current;
  double dt;
  time(&wtime_current);
  dt = difftime(wtime_current, _cs_timer_stdc_time_start);
  timer->wall_sec = floor(dt);
  timer->wall_nsec = (dt - timer->wall_sec) * 1.0e-9;
}

/*----------------------------------------------------------------------------
 * Null function to get CPU time.
 *
 * parameters:
 *   timer  <-- pointer to timing structure (ignored)
 *----------------------------------------------------------------------------*/

void
_cs_timer_cpu_null(cs_timer_t  *timer)
{
  CS_UNUSED(timer);
}

#if defined (HAVE_CLOCK_GETTIME)

/*----------------------------------------------------------------------------
 * Get CPU time using clock_gettime().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_cpu_clock_gettime(cs_timer_t  *timer)
{
  struct timespec cpu_time;
  (void)clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_time);
  timer->cpu_sec = cpu_time.tv_sec;
  timer->cpu_nsec = cpu_time.tv_nsec;
}

#endif /* defined (HAVE_CLOCK_GETTIME) */

#if defined (HAVE_GETRUSAGE)

/*----------------------------------------------------------------------------
 * Get CPU time using clock_getrusage().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_cpu_getrusage(cs_timer_t  *timer)
{
  struct rusage  usage;
  getrusage(RUSAGE_SELF, &usage);
  timer->cpu_sec = usage.ru_utime.tv_sec + usage.ru_stime.tv_sec;
  timer->cpu_nsec = (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec)*1000;
}

#endif /* defined (HAVE_GETRUSAGE) */

#if defined(_POSIX_SOURCE)

/*----------------------------------------------------------------------------
 * Get CPU time using times().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_cpu_times(cs_timer_t  *timer)
{
  struct tms  ptimer;
  clock_t ticks;
  times(&ptimer);
  ticks = ptimer.tms_utime + ptimer.tms_stime;
  timer->cpu_sec = ticks / _cs_timer_unit;
  timer->cpu_nsec = (double)(ticks % _cs_timer_unit)*1.e+9 / _cs_timer_unit;
 }

#endif /* defined(_POSIX_SOURCE) */

/*----------------------------------------------------------------------------
 * Get CPU time using clock().
 *
 * parameters:
 *   timer  <-- pointer to timing structure
 *----------------------------------------------------------------------------*/

void
_cs_timer_cpu_stdc_clock(cs_timer_t  *timer)
{
  clock_t clock_current = clock() - _cs_timer_clock_start;
  timer->cpu_sec = clock_current / CLOCKS_PER_SEC;
  timer->cpu_nsec = (clock_current % CLOCKS_PER_SEC)*1.e+9/CLOCKS_PER_SEC;
}

/*----------------------------------------------------------------------------
 * Initialize timers.
 *----------------------------------------------------------------------------*/

static void
_cs_timer_initialize(void)
{
  _cs_timer_start.wall_sec = 0;
  _cs_timer_start.wall_nsec = 0;
  _cs_timer_start.cpu_sec = 0;
  _cs_timer_start.cpu_nsec = 0;

  /* Select timing methods, trying highest resolution first */

#if defined (HAVE_CLOCK_GETTIME)

  struct timespec ts_time;

  if (_cs_timer_wall_method == CS_TIMER_DISABLE) {
    if (clock_gettime(CLOCK_REALTIME, &ts_time) == 0) {
      _cs_timer_start.wall_sec = ts_time.tv_sec;
      _cs_timer_start.wall_nsec = ts_time.tv_nsec;
      _cs_timer_wall_method = CS_TIMER_CLOCK_GETTIME;
      _cs_timer_wall = _cs_timer_wall_clock_gettime;
    }
  }

#if defined (HAVE_CLOCK_GETCPUCLOCKID)
  if (_cs_timer_cpu_method == CS_TIMER_DISABLE) {
    clockid_t clock_id;
    if (clock_getcpuclockid(0, &clock_id) == 0) {
      /* ENOENT could be allowed with process binding,
         but binding needs to be checked. */
      if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_time) == 0) {
        _cs_timer_cpu_method = CS_TIMER_CLOCK_GETTIME;
        _cs_timer_cpu = _cs_timer_cpu_clock_gettime;
      }
    }
  }
#endif

#endif

#if defined (HAVE_GETTIMEOFDAY)

  if (_cs_timer_wall_method == CS_TIMER_DISABLE) {
    static struct timeval  tv_time;
    if (gettimeofday(&tv_time, NULL) == 0) {
      _cs_timer_start.wall_sec = tv_time.tv_sec;
      _cs_timer_start.wall_nsec = tv_time.tv_usec*1000;
      _cs_timer_wall_method = CS_TIMER_GETTIMEOFDAY;
      _cs_timer_wall = _cs_timer_wall_gettimeofday;
    }
  }

#endif

#if defined(HAVE_GETRUSAGE)

  if (_cs_timer_cpu_method == CS_TIMER_DISABLE) {
    struct rusage  usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
      _cs_timer_cpu_method = CS_TIMER_GETRUSAGE;
      _cs_timer_cpu = _cs_timer_cpu_getrusage;
    }
  }

#endif

#if defined(_POSIX_SOURCE)

  if (_cs_timer_cpu_method == CS_TIMER_DISABLE) {
    static struct tms  ptimer;
    _cs_timer_unit = sysconf(_SC_CLK_TCK);
    if (_cs_timer_unit != -1 && times(&ptimer) != -1) {
      _cs_timer_cpu_method = CS_TIMER_TIMES;
      _cs_timer_cpu = _cs_timer_cpu_times;
    }
  }

#endif /* defined(_POSIX_SOURCE) */

  /* Use minimal C library functions */

  if (_cs_timer_wall_method == CS_TIMER_DISABLE) {
    time_t wtime_current;
    if (time(&wtime_current) != (time_t)-1) {
      _cs_timer_stdc_time_start = time(&wtime_current);
      _cs_timer_wall_method = CS_TIMER_TIME;
      _cs_timer_wall = _cs_timer_wall_stdc_time;
    }
  }

  if (_cs_timer_cpu_method == CS_TIMER_DISABLE) {
    _cs_timer_clock_start = clock();
    if (_cs_timer_clock_start != (clock_t)-1) {
      _cs_timer_cpu_method = CS_TIMER_CLOCK;
      _cs_timer_cpu = _cs_timer_cpu_stdc_clock;
    }
  }

  _cs_timer_initialized = true;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return Wall clock time
 *
 * \return elapsed time from first call of a function of the cs_timer_...()
 *         series, or -1 if unable to compute.
 */
/*----------------------------------------------------------------------------*/

double
cs_timer_wtime(void)
{
  cs_timer_t t1;

  /* Ensure initialization */

  if (_cs_timer_initialized == false)
    _cs_timer_initialize();

  /* Compute elapsed time */

  _cs_timer_wall(&t1);

  long long wall_nsec
    =  (t1.wall_sec - _cs_timer_start.wall_sec) * (long long)1000000000
      + t1.wall_nsec - _cs_timer_start.wall_nsec;

  return wall_nsec*1.e-9;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see cs_timer_cpu_time_method()), at least one of
 * the cs_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * \return current CPU time usage, or -1 if unable to compute.
 */
/*----------------------------------------------------------------------------*/

double
cs_timer_cpu_time(void)
{
  cs_timer_t t1;

  /* Ensure initialization */

  if (_cs_timer_initialized == 0)
    _cs_timer_initialize();

  /* Compute CPU time */

  _cs_timer_cpu(&t1);

  long long cpu_nsec
    =  (t1.cpu_sec - _cs_timer_start.cpu_sec) * (long long)1000000000
       + t1.cpu_nsec - _cs_timer_start.cpu_nsec;

  return cpu_nsec*1.e-9;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return separate user and system CPU times.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available, this function will return -1 values.
 *
 * \param [out] user_time current user CPU usage.
 * \param [out] system_time current system CPU usage.
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_cpu_times(double  *user_time,
                   double  *system_time)

{
  /* Ensure initialization */

  if (_cs_timer_initialized == 0)
     _cs_timer_initialize();

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

    if (_cs_timer_unit != -1 && times(&ptimer) != -1) {
      *user_time    = ((double)ptimer.tms_utime)  / _cs_timer_unit;
      *system_time  = ((double)ptimer.tms_stime)  / _cs_timer_unit;
    }

 }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a timer's value
 *
 * \return timer structure.
 */
/*----------------------------------------------------------------------------*/

cs_timer_t
cs_timer_time(void)
{
  cs_timer_t time_current;

  /* Ensure initialization */

  if (_cs_timer_initialized == false)
    _cs_timer_initialize();

  /* Compute elapsed time */

  _cs_timer_wall(&time_current);
  _cs_timer_cpu(&time_current);

  return time_current;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the difference between 2 timers.
 *
 * \param[in]  t0  oldest timer value
 * \param[in]  t1  most recent timer value
 *
 * \return last - first timer value.
 */
/*----------------------------------------------------------------------------*/

cs_timer_counter_t
cs_timer_diff(const cs_timer_t  *t0,
              const cs_timer_t  *t1)
{
  cs_timer_counter_t retval;

  retval.wall_nsec =  (t1->wall_sec - t0->wall_sec) * (long long)1000000000
                     + t1->wall_nsec - t0->wall_nsec;
  retval.cpu_nsec =  (t1->cpu_sec - t0->cpu_sec) * (long long)1000000000
                    + t1->cpu_nsec - t0->cpu_nsec;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return method used to return wall clock time.
 *
 * \return short description of method used to return wall clock time.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_timer_wtime_method(void)
{
  if (_cs_timer_initialized == 0)
    _cs_timer_initialize();

  switch(_cs_timer_wall_method) {

  case CS_TIMER_CLOCK_GETTIME:
    return _("clock_gettime() function");
  case CS_TIMER_GETTIMEOFDAY:
    return _("gettimeofday() function");
  case CS_TIMER_TIME:
    return _("Iso C time() function");
  default:
    return _("Disabled");

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return method used to return CPU time.
 *
 * \return short description of method used to return CPU time.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_timer_cpu_time_method(void)
{
  if (_cs_timer_initialized == 0)
    _cs_timer_initialize();

  switch(_cs_timer_wall_method) {

  case CS_TIMER_CLOCK_GETTIME:
    return _("clock_gettime() function");
  case CS_TIMER_GETRUSAGE:
    return _("getrusage() function");
  case CS_TIMER_TIMES:
    return _("Posix times() function");
  case CS_TIMER_CLOCK:
    return _("Iso C clock() function");
  case CS_TIMER_DISABLE:
    return _("Disabled");
  default:
    return _("Disabled");

  }
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
