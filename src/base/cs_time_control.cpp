/*============================================================================
 * Time dependency control for variables or properties.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_log.h"
#include "base/cs_map.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_time_control.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_control.cpp
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Default cs_time_control_f constructor
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t::cs_time_control_t()
{
  type = CS_TIME_CONTROL_TIME_STEP;
  at_start = false;
  at_first = false;
  at_end = false;
  start_nt = -1;
  end_nt = -1;
  interval_nt = 1;
  control_func = nullptr;
  control_input = nullptr;
  current_state = false;
  current_time_step = -1;
  last_nt = -2;
  last_t = -HUGE_VAL;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Time control constructor based on time step options.
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t::cs_time_control_t
(
  int    nt_start,     //<! start time step (or < 0 for unlimited)
  int    nt_end,       //<! end time step (or < 0 for unlimited)
  int    nt_interval,  //<! time step interval (< 0 if no periodic output)
  bool   at_start_,    //<! always active at start ?
  bool   at_end_,      //<! always active at end ?
  bool   at_first_     //<! always active at first step of this run ?
)
{
  *this = cs_time_control_t();

  at_start = at_start_;
  at_first = at_first_;
  at_end = at_end_;

  type = CS_TIME_CONTROL_TIME_STEP;

  if (nt_start < 0)
    nt_start = -1;
  if (nt_end < 0)
    nt_end = -1;
  if (nt_interval < 1)
    nt_interval = -1;

  start_nt = nt_start;
  end_nt = nt_end;
  interval_nt = nt_interval;
}

 /*----------------------------------------------------------------------------
 *!
 * \brief Time control constructor based on physical time options.
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t::cs_time_control_t
(
  double  t_start,     //<! start time step (or < 0 for unlimited)
  double  t_end,       //<! end time step (or < 0 for unlimited)
  double  t_interval,  //<! time interval (< 0 if no periodic output)
  bool    at_start_,   //<! always active at start ?
  bool    at_end_,     //<! always active at end ?
  bool    at_first_    //<! always active at first step of this run ?
)
{
  *this = cs_time_control_t();

  at_start = at_start_;
  at_first = at_first_;
  at_end = at_end_;

  type = CS_TIME_CONTROL_TIME;

  if (t_start < 0)
    t_start = -1;
  if (t_end < 0)
    t_end = -1;
  if (t_interval <= 0)
    t_interval = 0;

  start_t = t_start;
  end_t = t_end;
  interval_t = t_interval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Simple time control constructor based on external function.
 *
 * \remark: if the input pointer is non-null, it must point to valid data
 *          when the control function is called, so that value or structure
 *          should not be temporary (i.e. local);
 *
 * \param[in]  func        pointer to time control function.
 * \param[in]  input       pointer to optional (untyped) value or
 *                         structure, or nullptr.
 * \param[in]  at_start_   always active at start ?
 * \param[in]  at_end_     always active at end ?
 * \param[in]  at_first_   always active at first step of this run ?
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t::cs_time_control_t
(
  cs_time_control_func_t  *func,
  void                    *input,
  bool                     at_start_,
  bool                     at_end_,
  bool                     at_first_
)
{
  *this = cs_time_control_t();

  type = CS_TIME_CONTROL_FUNCTION;

  at_start = at_start_;
  at_first = at_first_;
  at_end = at_end_;

  start_nt = 0;
  end_nt = 0;
  interval_nt = 0;

  control_func = func;
  control_input = input;

  current_state = true;
  current_time_step = -1;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Indicate if a time control is active or not at the given time.
 *
 * If the time control or time step argument is nullptr, true is returned.
 *
 * \param[in]  ts  time step structure
 *
 * \return  true if active, false if inactive
 */
/*----------------------------------------------------------------------------*/

bool
cs_time_control_t::is_active(const cs_time_step_t  *ts)
{
  bool retval = false;

  if (ts == nullptr)
    retval = true;

  else {
    if (this->current_time_step == ts->nt_cur)
      retval = this->current_state;

    else {
      switch (this->type) {
      case CS_TIME_CONTROL_TIME_STEP:
        if (   this->interval_nt > 0
            && ts->nt_cur > ts->nt_prev
            && ts->nt_cur % (this->interval_nt) == 0)
          retval = true;
        if (this->start_nt > ts->nt_cur)
          retval = false;
        if (this->end_nt >= 0 && this->end_nt < ts->nt_cur)
          retval = false;
        break;

      case CS_TIME_CONTROL_TIME:
        if (this->interval_t > 0) {
          double delta_t = ts->t_cur - this->last_t;
          /* Ensure output is not spaced more than required interval
             (avoid missing output when time step varies) */
          if (delta_t >= this->interval_t * (1. + 1e-6))
            retval = true;

          /* Try to align output with frequency */
          else {
            double dt = ts->dt[0];
            double tp =   ts->t_cur
                        - this->interval_t*floor(ts->t_cur/this->interval_t);
            if (tp < dt && tp < (ts->t_cur - this->last_t))
              retval = true;
          }
          if (this->start_t > ts->t_cur)
            retval = false;
          if (this->end_t >= 0 && this->end_t < ts->nt_cur)
            retval = false;

        }
        break;

      case CS_TIME_CONTROL_FUNCTION:
        retval = this->control_func(ts, this->control_input);
      }

    }

    if (ts->nt_cur == ts->nt_prev && this->at_start)
      retval = true;
    if (ts->nt_cur == 1 && this->at_first)
      retval = true;
    if (ts->nt_cur == ts->nt_max && this->at_end)
      retval = true;

    if (this->current_time_step < ts->nt_cur) {
      this->current_time_step = ts->nt_cur;
      this->current_state = retval;
      if (retval) {
        this->last_nt = ts->nt_cur;
        this->last_t = ts->t_cur;
      }
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Get text description of time control configuration.
 *
 * If the time step argument is nullptr, true is returned.
 *
 * \param[out]  desc       description string
 * \param[in]   desc_size  description string maximum size
 *
 * \return  true if active, false if inactive
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_t::get_description
(
  char     *desc,
  size_t    desc_size
) const
{
  char b[256] = "";  /* should be more than enough */
  char *s = b;

  switch (this->type) {
  case CS_TIME_CONTROL_TIME_STEP:
    if (this->interval_nt == 1)
      s += sprintf(s, _(", every time step"));
    else if (this->interval_nt > 1)
      s += sprintf(s, _(", every %d time steps"), this->interval_nt);
    if (this->start_nt > 0)
      s += sprintf(s, _(", start %d"), this->start_nt);
    if (this->end_nt > 0)
      s += sprintf(s, _(", end %d"), this->end_nt);
    break;

  case CS_TIME_CONTROL_TIME:
    if (this->interval_t >= 0) {
      if (this->interval_t <= 0)
        s += sprintf(s, _(", every time step"));
      else
        s += sprintf(s, _(", every %g s"), this->interval_t);
    }
    if (this->start_t > 0)
      s += sprintf(s, _(", start %g s"), this->start_t);
    if (this->end_nt > 0)
      s += sprintf(s, _(", end %g s"), this->end_t);
    break;

  case CS_TIME_CONTROL_FUNCTION:
    s += sprintf(s, _(", function-based"));
  }

  if (this->at_start)
    s += sprintf(s, _(", at start"));
  if (this->at_end)
    s += sprintf(s, _(", at end"));

  int shift = 0;
  while (b[shift] == ' ' || b[shift] == ',')
    shift++;

  strncpy(desc, b+shift, desc_size);
  if (desc_size > 0)
    desc[desc_size-1] = '\0';
}

/*----------------------------------------------------------------------------
 *!
 * \brief Indicate if a time control is active or not at the given time.
 *
 * If the time control or time step argument is nullptr, true is returned.
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

  if (tc == nullptr)
    retval = true;
  else
    retval = tc->is_active(ts);

  return retval;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Simple time control initialization based on time step options.
 *
 * \param[in]  tc           pointer to time control structure.
 * \param[in]  nt_start     start time step (or < 0 for unlimited)
 * \param[in]  nt_end       end time step (or < 0 for unlimited)
 * \param[in]  nt_interval  time step interval (< 0 if no periodic output)
 * \param[in]  at_start     always active at start ?
 * \param[in]  at_start     always active at end ?
 */
/*----------------------------------------------------------------------------*/

[[deprecated("Use cs_time_control_t constructor instead")]] void
cs_time_control_init_by_time_step(cs_time_control_t  *tc,
                                  int                 nt_start,
                                  int                 nt_end,
                                  int                 nt_interval,
                                  bool                at_start,
                                  bool                at_end)
{
  *tc = cs_time_control_t(nt_start, nt_end, nt_interval,
                          at_start, at_end);
}

/*----------------------------------------------------------------------------
 *!
 * \brief Get text description of time control configuration.
 *
 * If the time control or time step argument is nullptr, true is returned.
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
  if (tc == nullptr)
    snprintf(desc, desc_size, "always active");
  else
    tc->get_description(desc, desc_size);

  if (desc_size > 0)
    desc[desc_size-1] = '\0';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy contents of a cs_time_control to another instance.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_copy
(
  const cs_time_control_t *src, /*!<[in]  Instance to copy data from */
  cs_time_control_t       *dst  /*!<[out] Instance to copy data to */
)
{
  dst->type = src->type;
  dst->at_start = src->at_start;
  dst->at_first = src->at_first;
  dst->at_end   = src->at_end;

  dst->start_nt = src->start_nt;
  dst->start_t  = src->start_t;

  dst->end_nt   = src->end_nt;
  dst->end_t    = src->end_t;

  dst->interval_nt = src->end_nt;
  dst->interval_t  = src->end_t;

  dst->control_func = src->control_func;
  dst->control_input = src->control_input;

  dst->current_state = src->current_state;
  dst->current_time_step = src->current_time_step;

  dst->last_nt = src->last_nt;
  dst->last_t = src->last_t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy contents of the default cs_time_control to another instance.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_control_copy_from_default
(
  cs_time_control_t *tc /*!<[out] Instance to copy data to */
)
{
  *tc = cs_time_control_t();
}

/*----------------------------------------------------------------------------*/
