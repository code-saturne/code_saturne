/*============================================================================
 * Time dependency control for variables or properties.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_map.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_control.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_control.c
        Time dependency control for variables or properties.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

static const cs_time_control_t  cs_time_control_default
= {
  .type = CS_TIME_CONTROL_TIME_STEP,
  .at_start = false,
  .at_end = false,
  .start_nt = -1,
  .end_nt = -1,
  .interval_nt = 1,
  .control_func = NULL,
  .control_input = NULL,
  .current_state = false,
  .current_time_step = -1,
  .last_nt = -2,
  .last_t = -HUGE_VAL
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Base initialization for time control.
 *
 * \param[in]  tc        pointer to time control structure.
 * \param[in]  at_start  always active at start ?
 * \param[in]  at_start  always active at end ?
 */
/*----------------------------------------------------------------------------*/

static void
_time_control_init_base(cs_time_control_t  *tc,
                        bool                at_start,
                        bool                at_end)
{
  memset(tc, 0, sizeof(cs_time_control_t));

  *tc = cs_time_control_default;

  tc->at_start = at_start;
  tc->at_end = at_end;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Indicate if a time control is active or not at the given time.
 *
 * If the time control or time step argument is NULL, true is returned.
 *
 * \param[in]  tc  time control structure
 * \param[in]  ts  time step structure
 *
 * \return  true if active, false if inactive
 */
/*----------------------------------------------------------------------------*/

bool
cs_time_control_is_active(cs_time_control_t     *tc,
                          const cs_time_step_t  *ts)
{
  bool retval = false;

  if (tc == NULL || ts == NULL)
    retval = true;

  else {
    if (tc->current_time_step == ts->nt_cur)
      retval = tc->current_state;

    else {
      switch (tc->type) {
      case CS_TIME_CONTROL_TIME_STEP:
        if (   tc->interval_nt > 0
            && ts->nt_cur > ts->nt_prev
            && ts->nt_cur % (tc->interval_nt) == 0)
          retval = true;
        if (tc->start_nt > ts->nt_cur)
          retval = false;
        if (tc->end_nt >= 0 && tc->end_nt < ts->nt_cur)
          retval = false;
        break;

      case CS_TIME_CONTROL_TIME:
        {
          double  delta_t = ts->t_cur - tc->last_t;
          if (   delta_t >= tc->interval_t*(1-1e-6)
              && tc->interval_t > 0)
            retval = true;
          if (tc->start_t > ts->t_cur)
            retval = false;
          if (tc->end_t >= 0 && tc->end_t < ts->nt_cur)
            retval = false;
        }
        break;

      case CS_TIME_CONTROL_FUNCTION:
        retval = tc->control_func(ts, tc->control_input);
      }

    }

    if (ts->nt_cur == ts->nt_prev && tc->at_start)
      retval = true;
    if (ts->nt_cur == ts->nt_max && tc->at_end)
      retval = true;

  }

  if (tc->current_time_step < ts->nt_cur) {
    tc->current_time_step = ts->nt_cur;
    tc->current_state = retval;
    if (retval) {
      tc->last_nt = ts->nt_cur;
      tc->last_t = ts->t_cur;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Simple time control initialization based on time step options.
 *
 * \param[in]  tc        pointer to time control structure.
 * \param[in]  nt_start  start time step (or < 0 for unlimited)
 * \param[in]  nt_end    end time step (or < 0 for unlimited)
 * \param[in]  at_start  always active at start ?
 * \param[in]  at_start  always active at end ?
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_init_by_time_step(cs_time_control_t  *tc,
                                  int                 nt_start,
                                  int                 nt_end,
                                  int                 nt_interval,
                                  bool                at_start,
                                  bool                at_end)
{
  _time_control_init_base(tc, at_start, at_end);

  tc->type = CS_TIME_CONTROL_TIME_STEP;

  if (nt_start < 0)
    nt_start = -1;
  if (nt_end < 0)
    nt_end = -1;
  if (nt_interval < 1)
    nt_interval = -1;

  tc->start_nt = nt_start;
  tc->end_nt = nt_end;
  tc->interval_nt = nt_interval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Simple time control initialization based on physical time options.
 *
 * \param[in]  tc        pointer to time control structure.
 * \param[in]  t_start   start time (or < 0 for unlimited)
 * \param[in]  t_end     end time (or < 0 for unlimited)
 * \param[in]  at_start  always active at start ?
 * \param[in]  at_start  always active at end ?
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_init_by_time(cs_time_control_t  *tc,
                             double              t_start,
                             double              t_end,
                             double              t_interval,
                             bool                at_start,
                             bool                at_end)
{
  _time_control_init_base(tc, at_start, at_end);

  tc->type = CS_TIME_CONTROL_TIME;

  if (t_start < 0)
    t_start = -1;
  if (t_end < 0)
    t_end = -1;
  if (t_interval <= 0)
    t_interval = 0;

  tc->start_t = t_start;
  tc->end_t = t_end;
  tc->interval_t = t_interval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Simple time control initialization based on external function.
 *
 * \remark: if the input pointer is non-NULL, it must point to valid data
 *          when the control function is called, so that value or structure
 *          should not be temporary (i.e. local);
 *
 * \param[in]  tc             pointer to time control structure.
 * \param[in]  control_func   pointer to time control funcction.
 * \param[in]  control_input  pointer to optional (untyped) value or structure,
 *                            or NULL.
 * \param[in]  at_start       always active at start ?
 * \param[in]  at_start       always active at end ?
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_init_by_func(cs_time_control_t       *tc,
                             cs_time_control_func_t  *control_func,
                             void                    *control_input,
                             bool                     at_start,
                             bool                     at_end)
{
  _time_control_init_base(tc, at_start, at_end);

  tc->type = CS_TIME_CONTROL_FUNCTION;

  tc->control_func = control_func;
  tc->control_input = control_input;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Get text description of time control configuration.
 *
 * If the time control or time step argument is NULL, true is returned.
 *
 * \param[in]   tc         time control structure
 * \param[out]  desc       description string
 * \param[in]   desc_size  description string maximum size
 *
 * \return  true if active, false if inactive
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_get_description(const cs_time_control_t  *tc,
                                char                     *desc,
                                size_t                    desc_size)
{
  char b[256] = "";  /* should be more than enough */
  char *s = b;

  if (tc == NULL) {
    snprintf(s, 256, "always active");
  }

  else {

    switch (tc->type) {
    case CS_TIME_CONTROL_TIME_STEP:
      if (tc->interval_nt == 1)
        s += sprintf(s, _(", every time step"));
      else if (tc->interval_nt > 1)
        s += sprintf(s, _(", every %d time steps"), tc->interval_nt);
      if (tc->start_nt > 0)
        s += sprintf(s, _(", start %d"), tc->start_nt);
      if (tc->end_nt > 0)
        s += sprintf(s, _(", end %d"), tc->end_nt);
      break;

      case CS_TIME_CONTROL_TIME:
        if (tc->interval_t >= 0) {
          if (tc->interval_t <= 0)
            s += sprintf(s, _(", every time step"));
          else
            s += sprintf(s, _(", every %g s"), tc->interval_t);
        }
        if (tc->start_t > 0)
          s += sprintf(s, _(", start %g s"), tc->start_t);
        if (tc->end_nt > 0)
          s += sprintf(s, _(", end %g s"), tc->end_t);
        break;

      case CS_TIME_CONTROL_FUNCTION:
        s += sprintf(s, _(", function-based"));
    }

    if (tc->at_start)
      s += sprintf(s, _(", at start"));
    if (tc->at_end)
      s += sprintf(s, _(", at end"));
  }

  int shift = 0;
  while (b[shift] == ' ' || b[shift] == ',')
    shift++;

  strncpy(desc, b+shift, desc_size);
  if (desc_size > 0)
    desc[desc_size-1] = '\0';
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
