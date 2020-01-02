#ifndef __CS_TIMER_STATS_H__
#define __CS_TIMER_STATS_H__

/*============================================================================
 * Application timer statistics and graphs
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "cs_timer.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize timer statistics handling.
 *
 * This creates 2 statistic timer trees, whose roots ids are:
 * - 0 for computational operations
 * - 1 for computational stages
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize timer statistics handling.
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a start time for time stats.
 *
 * This is useful to shift the time id for restarts. This function must
 * not be called after \ref cs_timer_stats_increment_time_step.
 *
 * \param[in]  time_id  associated starting time id
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_set_start_time(int time_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global timer statistics plot options.
 *
 * This function is only effective before the first call to
 * \ref cs_timer_stats_increment_time_step.
 *
 * \param[in]  format             associated file format
 * \param[in]  frequency          plot every n time steps
 * \param[in]  n_buffer_steps     number of time steps in output buffer if
 *                                file is not to be kept open
 * \param[in]  flush_wtime        elapsed time interval between file flushes
 *                                (if < 0, no forced flush)
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_set_plot_options(cs_time_plot_format_t   format,
                                int                     frequency,
                                int                     n_buffer_steps,
                                double                  flush_wtime);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Increment time step for timer statistics.
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_increment_time_step(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a timer statistics structure.
 *
 * \param[in]  parent_name  name of parent statistic, or NULL
 * \param[in]  name         associated canonical name
 * \param[in]  label        associated label, or NULL
 *
 * \return  id of new timer stats structure
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_create(const char  *parent_name,
                      const char  *name,
                      const char  *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of a defined statistic based on its name.
 *
 * If no timer with the given name exists, -1 is returned.
 *
 * \param[in]  name   statistic name
 *
 * \return  id of the statistic, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_id_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable or disable plotting for a timer statistic.
 *
 * By default plotting is enabled for all statistics, except root statistic 1
 * (as it measures the same total time as root 0, with a different subtree).
 *
 * This function is only effective before the first call to
 * \ref cs_timer_stats_increment_time_step.
 *
 * \param[in]  id    id of statistic
 * \param[in]  plot  0 to disable, 1 to enable
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_set_plot(int  id,
                        int  plot);

/*----------------------------------------------------------------------------*/
/*!
 * \brief indicate if a timer for a given statistic is currently active.
 *
 * \param[in]  id  id of statistic
 *
 * \return     1 if active, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_is_active(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start a timer for a given statistic.
 *
 * Parents of the current statistic are also started, if not active.
 *
 * If a timer with the same root but different parents is active, we assume
 * the current operation is a subset of the active timer, so the timer is
 * not started, so as to avoid having a sum of parts larger than the total.
 *
 * \param[in]  id  id of statistic
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_start(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Stop a timer for a given statistic.
 *
 * Children of the current statistic are also stopped, if active.
 *
 * \param[in]  id  id of statistic
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_stop(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start a timer for a given statistic, stopping previous timers
 *        of the same type which are not a parent, and starting inactive
 *        parent timers if necessary.
 *
 * \param[in]  id  id of statistic
 *
 * \return  id of previously active statistic, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_switch(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a timing range to an inactive timer.
 *
 * This does not modify parent timers, so consistency of active and inactive
 * timers must be ensured by the caller.
 *
 * \param[in]  id  id of statistic
 * \param[in]  t0  oldest timer value
 * \param[in]  t1  most recent timer value
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_add_diff(int  id,
                        const cs_timer_t    *t0,
                        const cs_timer_t    *t1);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define default timer statistics
 *
 * This creates 2 statistic timer trees, whose roots ids are:
 * - 0 for computational operations
 * - 1 for computational stages
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_define_defaults(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIMER_STATS_H__ */
