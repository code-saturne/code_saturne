#ifndef __CS_TIME_STEP_H__
#define __CS_TIME_STEP_H__

/*============================================================================
 * Base time step data.
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

  int           is_variable;  /* 0 if time step is fixed in time,
                                 1 if the time step is variable. */
  int           is_local;     /* 0 if time step is uniform in space,
                                 1 if it is local in space (in which case
                                 the time value is only a reference. */

  int           nt_prev;      /* absolute time step number reached by previous
                                 computation */
  int           nt_cur;       /* current absolute time step number */
  int           nt_max;       /* maximum absolute time step number */
  int           nt_ini;       /* Number of time steps for initialization */

  double        t_prev;       /* physical time reached by previous
                                 computation */
  double        t_cur;        /* current absolute time */
  double        t_max;        /* maximum absolute time */

  double        dt[3];        /* n, n-1, and n-2 time steps */
  double        dt_ref;       /* reference time step. */
  double        dt_next;      /* next (predicted) time step. */

} cs_time_step_t;


/* Time step options descriptor */
/*------------------------------*/

typedef struct {

  int       iptlro; /* Clip the time step with respect to the buoyant effects
                       - 0: false
                       - 1: true. */

  int       idtvar; /* Option for a variable time step
                       - -1: steady algorithm
                       -  0: constant time step
                       -  1: time step constant in space but variable in time
                       -  2: variable time step in space and in time. */

  double    coumax; /* Maximum Courant number (when idtvar is
                       different from 0). */

  double    cflmmx; /* Maximum Courant number for the continuity equation
                       in compressible model. */

  double    foumax; /* Maximum Fourier number (when idtvar is different from 0). */

  double    varrdt; /* Relative allowed variation of dt (when idtvar is
                       different from 0). */

  double    dtmin;  /* Minimum value of dt (when idtvar is different from 0).
                       Take
                       dtmin = min(ld/ud, sqrt(lt/(gdelta rho/rho)), ...). */

  double    dtmax;  /* Maximum value of dt (when idtvar is different from 0).
                       Take
                       dtmax = max(ld/ud, sqrt(lt/(gdelta rho/rho)), ...). */

  double    relxst; /* Relaxation coefficient for the steady algorithm. */

} cs_time_step_options_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main time step structure */

extern const cs_time_step_t  *cs_glob_time_step;

extern const cs_time_step_options_t  *cs_glob_time_step_options;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_time_step
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_time_step_t *
cs_get_glob_time_step(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_time_step_options
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_time_step_options_t *
cs_get_glob_time_step_options(void);

/*----------------------------------------------------------------------------
 * Define whether time step is variable or not
 *
 * parameters:
 *   is_variable <-- 0 if time step is variable in time, 1 if it is fixed
 *----------------------------------------------------------------------------*/

void
cs_time_step_define_variable(int  is_variable);

/*----------------------------------------------------------------------------
 * Define whether time step is local in space or not
 *
 * parameters:
 *   is_local <-- 0 if time step is uniform in space, 1 if it is local
 *----------------------------------------------------------------------------*/

void
cs_time_step_define_local(int  is_local);

/*----------------------------------------------------------------------------
 * Define maximum time step number
 *
 * parameters:
 *   nt_max <-- maximum time step number (unlimited if negative)
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

/*----------------------------------------------------------------------------*
 * Print the time stepping options to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_time_step_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TIME_STEP_H__ */
