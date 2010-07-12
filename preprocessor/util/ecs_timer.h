#ifndef __ECS_TIMER_H__
#define __ECS_TIMER_H__

/*============================================================================
 * Program timing information
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2009 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */

/* ECS library headers */

#include "cs_config.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Return Wall clock time
 *
 * returns:
 *   elapsed time from first call of a function of the ecs_timer_...()
 *   series, or -1 if unable to compute.
 */

double
ecs_timer_wtime(void);

/*
 * Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see ecs_timer_cpu_time_method()), at least one of
 * the ecs_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * returns:
 *   current CPU time usage, or -1 if unable to compute.
 */

double
ecs_timer_cpu_time(void);

/*
 * Return separate user and system CPU times.
 *
 * parameters:
 *   user_time   --> current user CPU usage.
 *   system_time --> current system CPU usage.
 */

void
ecs_timer_cpu_times(double *user_time,
                    double *system_time);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __ECS_TIMER_H__ */
