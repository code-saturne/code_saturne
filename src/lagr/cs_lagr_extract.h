#ifndef __CS_LAGR_EXTRACT_H__
#define __CS_LAGR_EXTRACT_H__

/*============================================================================
 * Extract information from lagrangian particles.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_lagr_particle.h"
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
 * Get the local number of particles.
 *
 * returns:
 *   current number of particles.
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_get_n_particles(void);

/*----------------------------------------------------------------------------
 * Extract a list of particles using an optional cell filter and
 * statistical density filter.
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * parameters:
 *   n_cells         <-- number of cells in filter
 *   cell_list       <-- optional list of containing cells filter
 *                       (1 to n numbering)
 *   density         <-- if < 1, fraction of particles to select
 *   n_particles     --> number of selected particles, or NULL
 *   particle_list   --> particle_list (1 to n numbering), or NULL
 *----------------------------------------------------------------------------*/

void
cs_lagr_get_particle_list(cs_lnum_t         n_cells,
                          const cs_lnum_t   cell_list[],
                          double            density,
                          cs_lnum_t        *n_particles,
                          cs_lnum_t        *particle_list);

/*----------------------------------------------------------------------------
 * Extract values for a set of particles.
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * parameters:
 *   particle_set  <-- associated particle set
 *   attr          <-- attribute whose values are required
 *   datatype      <-- associated value type
 *   stride        <-- number of values per particle
 *   component_id  <-- if -1 : extract the whole attribute
 *                     if >0 : id of the component to extract
 *   n_particles   <-- number of particles in filter
 *   particle_list <-- particle_list (1 to n numbering), or NULL
 *   values        --> particle values for given attribute
 *
 * returns:
 *   0 in case of success, 1 if attribute is not present
 *----------------------------------------------------------------------------*/

int
cs_lagr_get_particle_values(const cs_lagr_particle_set_t  *particles,
                            cs_lagr_attribute_t            attr,
                            cs_datatype_t                  datatype,
                            int                            stride,
                            int                            component_id,
                            cs_lnum_t                      n_particles,
                            const cs_lnum_t                particle_list[],
                            void                          *values);

/*----------------------------------------------------------------------------
 * Extract trajectory values for a set of particles.
 *
 * Trajectories are defined as a mesh of segments, whose start and end
 * points are copied in an interleaved manner in the segment_values array
 * (p1_old, p1_new, p2_old, p2_new, ... pn_old, pn_new).
 *
 * The output array must have been allocated by the caller and be of
 * sufficient size.
 *
 * parameters:
 *   particles      <-- associated particle set
 *   attr           <-- attribute whose values are required
 *   datatype       <-- associated value type
 *   stride         <-- number of values per particle
 *   component_id   <-- if -1 : extract the whole attribute
 *                      if >0 : id of the component to extract
 *   n_particles    <-- number of particles in filter
 *   particle_list  <-- particle_list (1 to n numbering), or NULL
 *   segment_values --> particle segment values
 *
 * returns:
 *    0 in case of success, 1 if attribute is not present
 *----------------------------------------------------------------------------*/

int
cs_lagr_get_trajectory_values(const cs_lagr_particle_set_t  *particles,
                              cs_lagr_attribute_t            attr,
                              cs_datatype_t                  datatype,
                              int                            stride,
                              int                            component_id,
                              cs_lnum_t                      n_particles,
                              const cs_lnum_t                particle_list[],
                              void                          *segment_values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_EXTRACT_H__ */
