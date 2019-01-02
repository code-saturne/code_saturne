#ifndef __CS_TIMER_H__
#define __CS_TIMER_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/* Information structure for precise timings */

typedef struct {

  long long    wall_sec;       /* wall-time seconds */
  long long    wall_nsec;      /* wall-time nanoseconds */
  long long    cpu_sec;        /* CPU time seconds */
  long long    cpu_nsec;       /* CPU time nanoseconds */

} cs_timer_t;

/* Information structure for timing counters */

typedef struct {

  long long    wall_nsec;      /* wall-time nanoseconds */
  long long    cpu_nsec;       /* CPU time nanoseconds */

} cs_timer_counter_t;

/*============================================================================
 * Public macros
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize timer counter.
 *
 * parameters:
 *   _t --> resulting counter.
 *----------------------------------------------------------------------------*/

#define CS_TIMER_COUNTER_INIT(_t)       \
  (_t.wall_nsec = 0,  \
   _t.cpu_nsec  = 0)

/*----------------------------------------------------------------------------
 * Add timer counter.
 *
 * The result may be identical to one of the 2 counters to add.
 *
 * parameters:
 *   _res --> resulting counter.
 *   _c0  <-- counter to add.
 *   _c1  <-- counter to add.
 *----------------------------------------------------------------------------*/

#define CS_TIMER_COUNTER_ADD(_res, _c0, _c1)       \
  (_res.wall_nsec = _c0.wall_nsec + _c1.wall_nsec,  \
   _res.cpu_nsec  = _c0.cpu_nsec + _c1.cpu_nsec)

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return Wall clock time
 *
 * returns:
 *   elapsed time from first call of a function of the cs_timer_...()
 *   series, or -1 if unable to compute.
 *----------------------------------------------------------------------------*/

double
cs_timer_wtime(void);

/*----------------------------------------------------------------------------
 * Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see cs_timer_cpu_time_method()), at least one of
 * the cs_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * returns:
 *   current CPU time usage, or -1 if unable to compute.
 *----------------------------------------------------------------------------*/

double
cs_timer_cpu_time(void);

/*----------------------------------------------------------------------------
 * Return separate user and system CPU times.
 *
 * parameters:
 *   user_time   --> current user CPU usage.
 *   system_time --> current system CPU usage.
 *----------------------------------------------------------------------------*/

void
cs_timer_cpu_times(double  *user_time,
                   double  *system_time);

/*----------------------------------------------------------------------------
 * Return a timer's value
 *
 * returns:
 *   timer structure.
 *----------------------------------------------------------------------------*/

cs_timer_t
cs_timer_time(void);

/*----------------------------------------------------------------------------
 * Compute the difference between 2 timers.
 *
 * parameters:
 *   t0 <-- oldest timer value
 *   t1 <-- most recent timer value
 *
 * returns:
 *   last - first timer value.
 *----------------------------------------------------------------------------*/

cs_timer_counter_t
cs_timer_diff(const cs_timer_t  *t0,
              const cs_timer_t  *t1);

/*----------------------------------------------------------------------------
 * Add the the difference between 2 timers to a counter.
 *
 * parameters:
 *   tc <-> pointer to timer counter
 *   t0 <-- oldest timer value
 *   t1 <-- most recent timer value
 *
 * returns:
 *   last - first timer value.
 *----------------------------------------------------------------------------*/

void
cs_timer_counter_add_diff(cs_timer_counter_t  *tc,
                          const cs_timer_t    *t0,
                          const cs_timer_t    *t1);

/*----------------------------------------------------------------------------
 * Return method used to return wall clock time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available, this function will return -1 values.
 *
 * returns:
 *   short description of method used to return wall clock time.
 *----------------------------------------------------------------------------*/

const char *
cs_timer_wtime_method(void);

/*----------------------------------------------------------------------------
 * Return method used to return CPU time.
 *
 * returns:
 *   short description of method used to return CPU time.
 *----------------------------------------------------------------------------*/

const char *
cs_timer_cpu_time_method(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIMER_H__ */
