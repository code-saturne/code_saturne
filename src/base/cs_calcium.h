#ifndef __CS_CALCIUM_H__
#define __CS_CALCIUM_H__

/*============================================================================
 * Basic CALCIUM-mappable functions for code coupling
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Type Definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read values, blocking until they are available.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_read_int(int                    rank_id,
                    int                   *iteration,
                    const char            *var_name,
                    int                    n_val_max,
                    int                   *n_val_read,
                    int                    val[]);

int
cs_calcium_read_double(int                    rank_id,
                       int                   *iteration,
                       const char            *var_name,
                       int                    n_val_max,
                       int                   *n_val_read,
                       double                 val[]);

/*----------------------------------------------------------------------------
 * Write values.
 *
 * parameters:
 *   rank_id    <-- communicating MPI rank id
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_int(int                    rank_id,
                     int                    iteration,
                     const char            *var_name,
                     int                    n_val,
                     const int              val[]);

int
cs_calcium_write_double(int                    rank_id,
                        int                    iteration,
                        const char            *var_name,
                        int                    n_val,
                        const double           val[]);

/*----------------------------------------------------------------------------
 * Set the CALCIUM-mappable function's verbosity
 *
 * parameters:
 *   n_echo <-- verbosity (none if -1, headers if 0,
 *              headers + n first and last elements if > 0.
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_verbosity(int  n_echo);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CALCIUM_H__ */

