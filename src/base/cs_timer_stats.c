/*============================================================================
 * Application timer statistics and graphs
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
#include <string.h>
#include <time.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_map.h"
#include "cs_timer.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_timer_stats.c
        Application timer statistics and graphs.

  These statistics are intended to provide a synthetic view of the relative
  costs of various operations and stages, and are defined as sets of trees
  so as to allow a form of grouping.

  This logic does not replace having specific timers for operators,
  which allow for more detail, but more difficult to provide an overview for.
  Timer statistics also allow for incrementing results from base timers
  (in addition to starting/stopping their own timers), so they may be used
  to assist logging and plotting of other timers.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/* Field key definitions */

typedef struct {

  char                *label;           /* Associated label */

  int                  root_id;         /* Parent id */
  int                  parent_id;       /* Parent id */

  bool                 plot;            /* true if plot desired */
  bool                 active;          /* true if active, false otherwise */
  cs_timer_t           t_start;         /* Start time if active */

  cs_timer_counter_t   t_cur;           /* Counter since last output */
  cs_timer_counter_t   t_tot;           /* Total time counter */

} cs_timer_stats_t;

/*-------------------------------------------------------------------------------
 * Local macro documentation
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

/* timer options */

static int                    _plot_frequency = 1;
static int                    _plot_buffer_steps = -1;
static double                 _plot_flush_wtime = 3600;
static cs_time_plot_format_t  _plot_format = CS_TIME_PLOT_CSV;

/* Timer status */

static int  _time_id = -1, _start_time_id = -1;
static int  _n_roots = 0;
static int  *_active_id = NULL;
static cs_time_plot_t  *_time_plot = NULL;

/* Field definitions */

static int  _n_stats = 0;
static int  _n_stats_max = 0;
static cs_timer_stats_t  *_stats= NULL;

static cs_map_name_to_id_t  *_name_map = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if a timer is a parent of another
 *
 * parameters:
 *   id_0  <-- id of statistic 0
 *   id_1  <-- id of statistic 1
 *
 * return:
 *   true if id_0 is a parent of id_1
 *----------------------------------------------------------------------------*/

static inline bool
_is_parent(int  id_0,
           int  id_1)
{
  bool retval = false;

  if (id_0 < 0 || id_0 == id_1)
    retval = true;
  else if (id_0 > id_1)
    retval = false;

  else {
    cs_timer_stats_t  *s = _stats + id_1;
    while (s->parent_id > -1) {
      if (s->parent_id == id_0) {
        retval = true;
        break;
      }
      else
        s = _stats + s->parent_id;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return common parent id of two statistics
 *
 * parameters:
 *   id_0  <-- id of statistic 0
 *   id_1  <-- id of statistic 1
 *
 * return:
 *   common parent id of two statistics
 *----------------------------------------------------------------------------*/

static inline int
_common_parent_id(int  id_0,
                  int  id_1)
{
  int p0 = id_0, p1 = id_1;

  while (p0 != p1 && p0 > -1 && p1 > -1) {
    if (p0 < p1)
      p1 = (_stats + p1)->parent_id;
    else
      p0 = (_stats + p0)->parent_id;
  }

  if (p0 != p1)
    p0 = -1;

  return p0;
}

/*----------------------------------------------------------------------------
 * Create time plots
 *----------------------------------------------------------------------------*/

static void
_build_time_plot(void)
{
  const char **stats_labels;
  BFT_MALLOC(stats_labels, _n_stats, const char *);

  int stats_count = 0;

  for (int stats_id = 0; stats_id < _n_stats; stats_id++) {
    cs_timer_stats_t  *s = _stats + stats_id;
    if (s->plot) {
      stats_labels[stats_count] = s->label;
      stats_count++;
    }
  }

  if (stats_count > 0)
    _time_plot = cs_time_plot_init_probe("timer_stats",
                                         "",
                                         _plot_format,
                                         true,
                                         _plot_flush_wtime,
                                         _plot_buffer_steps,
                                         stats_count,
                                         NULL,
                                         NULL,
                                         stats_labels);

  BFT_FREE(stats_labels);
}

/*----------------------------------------------------------------------------
 * Output time plots
 *----------------------------------------------------------------------------*/

static void
_output_time_plot(void)
{
  cs_real_t *vals;
  BFT_MALLOC(vals, _n_stats, cs_real_t);

  int stats_count = 0;

  for (int stats_id = 0; stats_id < _n_stats; stats_id++) {

    cs_timer_stats_t  *s = _stats + stats_id;
    if (s->plot) {
      vals[stats_count] = s->t_cur.wall_nsec*1e-9;
      stats_count++;
    }

  }

  cs_time_plot_vals_write(_time_plot,
                          _time_id,
                          -1.,
                          stats_count,
                          vals);

  BFT_FREE(vals);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
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
cs_timer_stats_initialize(void)
{
  int id;

  _time_id = 0;
  _start_time_id = 0;

  if (_name_map != NULL)
    cs_timer_stats_finalize();

  _name_map = cs_map_name_to_id_create();

  id = cs_timer_stats_create(NULL, "operations", "total");
  cs_timer_stats_start(id);

  id = cs_timer_stats_create(NULL, "stages", "total");
  cs_timer_stats_start(id);
  cs_timer_stats_set_plot(id, 0);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize timer statistics handling.
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_finalize(void)
{
  cs_timer_stats_increment_time_step();

  if (_time_plot != NULL)
    cs_time_plot_finalize(&_time_plot);

  _time_id = -1;

  for (int stats_id = 0; stats_id < _n_stats; stats_id++) {
    cs_timer_stats_t  *s = _stats + stats_id;
    BFT_FREE(s->label);
  }

  BFT_FREE(_stats);

  BFT_FREE(_active_id);
  _n_roots = 0;

  cs_map_name_to_id_destroy(&_name_map);

  _n_stats = 0;
  _n_stats_max = 0;
}

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
cs_timer_stats_set_start_time(int time_id)
{
  int id;

  if (_time_id <= 0 && _start_time_id <= 0) {
    _time_id = time_id;
    _start_time_id = time_id;
  }
}

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
                                double                  flush_wtime)
{
  _plot_format = format;

  _plot_frequency =  frequency;
  _plot_buffer_steps = n_buffer_steps;
  _plot_flush_wtime = flush_wtime;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Increment time step for timer statistics.
 */
/*----------------------------------------------------------------------------*/

void
cs_timer_stats_increment_time_step(void)
{
  cs_timer_t t_incr = cs_timer_time();

  /* Update start and current time for active statistics
     (should be only root statistics if used properly) */

  for (int stats_id = 0; stats_id < _n_stats; stats_id++) {
    cs_timer_stats_t  *s = _stats + stats_id;
    if (s->active) {
      cs_timer_counter_add_diff(&(s->t_cur), &(s->t_start), &t_incr);
      s->t_start = t_incr;
    }
  }

  /* Now output data */

  if (   _time_plot == NULL && _time_id < _start_time_id + 1
      && cs_glob_rank_id < 1)
    _build_time_plot();

  if (_time_id % _plot_frequency == 0) {

    if (_time_plot != NULL)
      _output_time_plot();

    for (int stats_id = 0; stats_id < _n_stats; stats_id++) {
      cs_timer_stats_t  *s = _stats + stats_id;
      CS_TIMER_COUNTER_ADD(s->t_tot, s->t_tot, s->t_cur);
      CS_TIMER_COUNTER_INIT(s->t_cur);
    }

  }

  _time_id += 1;
}

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
                      const char  *label)
{
  /* Determine parent id, create new series if none */

  const char *_parent_name = NULL;
  int parent_id = -1;
  int root_id = -1;

  if (parent_name != NULL) {
    if (strlen(parent_name) > 0)
      _parent_name = parent_name;
  }

  if (_parent_name == NULL) {
    BFT_REALLOC(_active_id, _n_roots+1, int);
    _active_id[_n_roots] = -1;
    root_id = _n_roots;
    _n_roots += 1;
  }
  else {
    parent_id = cs_map_name_to_id_try(_name_map, parent_name);
    if (parent_id < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Timer statistics \"%s\"\n"
                  " parent \"%s\" not defined."),
                name, parent_name);
  }

  cs_timer_stats_t  *s;

  /* Insert entry in map */

  int stats_id = cs_map_name_to_id(_name_map, name);

  if (stats_id < _n_stats) {
    s =  _stats + stats_id;
    bft_error(__FILE__, __LINE__, 0,
              _("Timer statistics \"%s\"\n"
                " is already defined, with id %d and parent %d."),
              name, stats_id, s->parent_id);
  }
  else
    _n_stats = stats_id + 1;

  /* Reallocate pointers if necessary */

  if (_n_stats > _n_stats_max) {
    if (_n_stats_max == 0)
      _n_stats_max = 8;
    else
      _n_stats_max *= 2;
    BFT_REALLOC(_stats, _n_stats_max, cs_timer_stats_t);
  }

  /* Now build new statistics */

  s =  _stats + stats_id;

  s->label = NULL;
  if (label != NULL) {
    size_t l_len = strlen(label);
    if (l_len > 0) {
      BFT_MALLOC(s->label, l_len + 1, char);
      strcpy(s->label, label);
    }
  }
  if (s->label == NULL) {
    BFT_MALLOC(s->label, strlen(name) + 1, char);
    strcpy(s->label, name);
  }

  s->parent_id = parent_id;
  if (root_id < 0)
    s->root_id = (_stats + parent_id)->root_id;
  else
    s->root_id = root_id;

  s->plot = true;
  s->active = false;

  CS_TIMER_COUNTER_INIT(s->t_cur);
  CS_TIMER_COUNTER_INIT(s->t_tot);

  return stats_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of a defined statistic based on its name.
 *
 * If no timer with the given name exists, -1 is returned.
 *
 * \param[in]  name   statitic name
 *
 * \return  id of the statistic, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_id_by_name(const char  *name)
{
  return cs_map_name_to_id_try(_name_map, name);
}

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
                        int  plot)
{
  if (id < 0 || id > _n_stats || _time_plot != NULL) return;

  cs_timer_stats_t  *s = _stats + id;
  s->plot = (plot != 0) ? true : false;
}

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
cs_timer_stats_is_active(int  id)
{
  int retval = 0;
  if (id >= 0 && id < _n_stats) {
    cs_timer_stats_t  *s = _stats + id;
    if (s->active)
      retval = 1;
  }
  return retval;
}

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
cs_timer_stats_start(int  id)
{
  if (id < 0 || id > _n_stats) return;

  cs_timer_stats_t  *s = _stats + id;

  cs_timer_t t_start = cs_timer_time();

  const int root_id = s->root_id;

  /* If a timer with the same root but different parents is active,
     simply return */

  if (! _is_parent(_active_id[root_id], id))
    return;

  int parent_id = _common_parent_id(id, _active_id[root_id]);

  /* Start timer and inactive parents */

  for (int p_id = id; p_id > parent_id; p_id = (_stats + p_id)->parent_id) {

    s = _stats + p_id;

    if (s->active == false) {
      s->active = true;
      s->t_start = t_start;
    }

  }

  _active_id[root_id] = id;
}

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
cs_timer_stats_stop(int  id)
{
  if (id < 0 || id > _n_stats) return;

  cs_timer_stats_t  *s = _stats + id;

  cs_timer_t t_stop = cs_timer_time();

  /* Stop timer and active children */

  const int root_id = s->root_id;

  while (_is_parent(id, _active_id[root_id])) {

    s = _stats + _active_id[root_id];

    if (s->active == true) {
      s->active = false;
      _active_id[root_id] = s->parent_id;
      cs_timer_counter_add_diff(&(s->t_cur), &(s->t_start), &t_stop);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start a timer for a given statistic, stopping previous timers
 *        of the same type which are not a parent, and starting inactive
 *        parent timers if necessary.
 *
 * \param[in]  id  id of statistic with same root
 *
 * \return  id of previously active statistic, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_timer_stats_switch(int  id)
{
  int retval = -1;

  if (id < 0 || id > _n_stats) return retval;

  cs_timer_stats_t  *s = _stats + id;

  cs_timer_t t_switch = cs_timer_time();

  const int root_id = s->root_id;

  retval = _active_id[root_id];

  if (_active_id[root_id] == id)
    return retval; /* Nothing to do, already current */

  int parent_id = _common_parent_id(id, _active_id[root_id]);

  /* Stop all active timers of same type which are lower level than the
     common parent. */

  while (parent_id != _active_id[root_id]) {

    s = _stats + _active_id[root_id];

    if (s->active == true) {
      s->active = false;
      _active_id[root_id] = s->parent_id;
      cs_timer_counter_add_diff(&(s->t_cur), &(s->t_start), &t_switch);
    }

  }

  /* Start all inactive timers of the same type which are lower level
     than the common parent */

  for (int p_id = id; p_id > parent_id; p_id = (_stats + id)->parent_id) {

    s = _stats + p_id;

    if (s->active == false) {
      s->active = true;
      s->t_start = t_switch;
    }

  }

  _active_id[root_id] = id;

  return retval;
}

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
cs_timer_stats_add_diff(int                id,
                        const cs_timer_t  *t0,
                        const cs_timer_t  *t1)
{
  if (id < 0 || id > _n_stats) return;

  cs_timer_stats_t  *s = _stats + id;

  if (s->active == false)
    cs_timer_counter_add_diff(&(s->t_cur), t0, t1);
}

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
cs_timer_stats_define_defaults(void)
{
  int id;

  /* Operations */

  cs_timer_stats_create("operations",
                        "mesh_processing",
                        "mesh processing");

  id = cs_timer_stats_create("mesh_processing",
                             "mesh_io",
                             "mesh io");
  cs_timer_stats_set_plot(id, 0);

  id = cs_timer_stats_create("operations",
                             "postprocessing_output",
                             "post-processing output");
  cs_timer_stats_set_plot(id, 0);

  /* Stages */

  id = cs_timer_stats_create ("stages",
                              "checkpoint_restart_stage",
                              "checkpoint/restart");

  id = cs_timer_stats_create("stages",
                             "postprocessing_stage",
                             "post-processing");
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
