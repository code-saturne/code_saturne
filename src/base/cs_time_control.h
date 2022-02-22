#ifndef __CS_TIME_CONTROL_H__
#define __CS_TIME_CONTROL_H__

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function pointer to a time control function.
 *
 * Each function of this sort may be used to determine whether a given
 * set of operations (such as variable or property updates) should
 * be done at the current time step.
 *
 * \remark: if the input pointer is non-NULL, it must point to valid data
 *          when the control function is called, so that value or structure
 *          should not be temporary (i.e. local);
 *
 * \param[in]  ts     current time step structure
 * \param[in]  input  pointer to optional (untyped) value or structure, or NULL.
 *
 * \return  true if active at given time step, false oltherwise
 */
/*----------------------------------------------------------------------------*/

typedef bool
(cs_time_control_func_t) (const cs_time_step_t  *ts,
                          void                  *input);

/* Time control types */
/*--------------------*/

/* Datatype enumeration */

typedef enum {

  CS_TIME_CONTROL_TIME_STEP,     /*!< control based on time step */
  CS_TIME_CONTROL_TIME,          /*!< control based on simulated time */
  CS_TIME_CONTROL_FUNCTION,      /*!< control based on function */

} cs_time_control_type_t;

/*----------------------------------------------------------------------------
 * Time control structure
 *----------------------------------------------------------------------------*/

typedef struct {

  /* Type and control parameters */

  cs_time_control_type_t  type;   /* control type */

  bool       at_start;            /* always active at start ? */
  bool       at_end;              /* always active at end ? */

  union {
    int      start_nt;            /* update start time step */
    double   start_t;             /* update start physical time */
  };

  union {
    int      end_nt;              /* update end time step */
    double   end_t;               /* update end physical time */
  };

  union {
    int      interval_nt;         /* interval, in time steps */
    double   interval_t;          /* interval, in physical time units */
  };

  cs_time_control_func_t  *control_func;   /* function for advanced control */
  void                    *control_input;  /* input associated function */

  /* Current state */

  bool    current_state;          /* query return value of current time step */
  int     current_time_step;      /* last time step queried */

  int     last_nt;                /* last active time step */
  double  last_t;                 /* last active physical time */

}  cs_time_control_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
                          const cs_time_step_t  *ts);

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
                                  bool                at_end);

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
                             bool                at_end);

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
                             bool                     at_end);

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
                                size_t                    desc_size);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_CONTROL_H__ */
