#ifndef __CS_TIME_STEP_H__
#define __CS_TIME_STEP_H__

/*============================================================================
 * Base time step data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* time step descriptor */
/*----------------------*/

typedef struct {

  int           nt_prev;      /* absolute time step number reached by previous
                                 computation */
  int           nt_cur;       /* current absolute time step number */
  int           nt_max;       /* maximum absolute time step number */

  double        t_prev;       /* physical time reached by previous
                                 computation */
  double        t_cur;        /* current absolute time */
  double        t_max;        /* maximum absolute time */

} cs_time_step_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main time step structure */

extern const cs_time_step_t  *cs_glob_time_step;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define maximum time step number
 *
 * parameters:
 *   nt_max  maximum time step number (unlimited if negative)
 *----------------------------------------------------------------------------*/

void
cs_time_step_define_nt_max(int  nt_max);

/*----------------------------------------------------------------------------
 * Define maximum time value
 *
 * parameters:
 *   t_max <-- maximum time value (unlimited if negative)
 *----------------------------------------------------------------------------*/

void
cs_time_step_define_t_max(double  t_max);

/*----------------------------------------------------------------------------
 * Set time values from previous (usually restarted) calculations
 *
 * parameters:
 *   nt_prev <-- previous time step number
 *   t_prev  <-- previous physical time
 *----------------------------------------------------------------------------*/

void
cs_time_step_define_prev(int     nt_prev,
                         double  t_prev);

/*----------------------------------------------------------------------------
 * Increment the global time step.
 *
 * parameters:
 *   dt <-- time step value to increment
 *----------------------------------------------------------------------------*/

void
cs_time_step_increment(double  dt);

/*----------------------------------------------------------------------------
 * Redefine the current time values.
 *
 * Remark: Using cs_time_step_increment() is preferred, but this function
 *         may be required for reverting to a previous time step.
 *
 * parameters:
 *   nt_cur <-- current time step number
 *   t_cur  <-- current physical time
 *----------------------------------------------------------------------------*/

void
cs_time_step_redefine_cur(int     nt_cur,
                          double  t_cur);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_STEP_H__ */
