#ifndef __CS_LOG_ITERATION_H__
#define __CS_LOG_ITERATION_H__

/*============================================================================
 * Log field and other array statistics at relevant time steps.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Add or update array not saved as permanent field to iteration log.
 *
 * parameters:
 *   name         <-- array name
 *   category     <-- category name
 *   loc_id       <-- associated mesh location id
 *   is_intensive <-- are the matching values intensive ?
 *   dim          <-- associated dimension (interleaved)
 *   val          <-- associated values
 *----------------------------------------------------------------------------*/

void
cs_log_iteration_add_array(const char                     *name,
                           const char                     *category,
                           const cs_mesh_location_type_t   loc_id,
                           bool                            is_intensive,
                           int                             dim,
                           const cs_real_t                 val[]);

/*----------------------------------------------------------------------------
 * Add or update clipping info for a given array.
 *
 * parameters:
 *   name         <-- array name
 *   dim          <-- associated dimension
 *   n_clip_min   <-- number of local clippings to minimum value
 *   n_clip_max   <-- number of local clippings to maximum value
 *   min_pre_clip <-- minimum values prior to clipping
 *   max_pre_clip <-- maximum values prior to clipping
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_clipping(const char       *name,
                          int               dim,
                          cs_lnum_t         n_clip_min,
                          cs_lnum_t         n_clip_max,
                          const cs_real_t   min_pre_clip[],
                          const cs_real_t   max_pre_clip[]);

/*----------------------------------------------------------------------------
 * Add or update clipping info for a field.
 *
 * parameters:
 *   f_id         <-- associated field id
 *   n_clip_min   <-- number of local clippings to minimum value
 *   n_clip_max   <-- number of local clippings to maximum value
 *   min_pre_clip <-- minimum values prior to clipping
 *   max_pre_clip <-- maximum values prior to clipping
 */
/*----------------------------------------------------------------------------*/

void
cs_log_iteration_clipping_field(int               f_id,
                                cs_lnum_t         n_clip_min,
                                cs_lnum_t         n_clip_max,
                                const cs_real_t   min_pre_clip[],
                                const cs_real_t   max_pre_clip[],
                                cs_lnum_t         n_clip_min_comp[],
                                cs_lnum_t         n_clip_max_comp[]);

/*----------------------------------------------------------------------------
 * Log L2 time residual for every variable fields.
 *----------------------------------------------------------------------------*/

void
cs_log_l2residual(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LOG_ITERATION_H__ */
