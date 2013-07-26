#ifndef __CS_LOG_ITERATION_H__
#define __CS_LOG_ITERATION_H__

/*============================================================================
 * Log field and other array statistics at relevant time steps.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Free arrays possible used by logging of array statistics.
 *----------------------------------------------------------------------------*/

void
cs_log_iteration_destroy_all(void);

/*----------------------------------------------------------------------------
 * Log field and other array statistics for the current time step.
 *----------------------------------------------------------------------------*/

void
cs_log_iteration(void);

/*----------------------------------------------------------------------------
 * Add array not saved as permanent field to logging of fields.
 *
 * parameters:
 *   name         <-- array name
 *   category     <-- category name
 *   loc_id       <-- associated mesh location id
 *   is_intensive <-- are the matching values intensive ?
 *   dimension    <-- associated dimension (interleaved)
 *   val          <-- associated values
 *----------------------------------------------------------------------------*/

void
cs_log_iteration_add_array(const char                     *name,
                           const char                     *category,
                           const cs_mesh_location_type_t   loc_id,
                           bool                            is_intensive,
                           int                             dim,
                           const cs_real_t                 val[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LOG_ITERATION_H__ */
