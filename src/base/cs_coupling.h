#ifndef __CS_COUPLING_H__
#define __CS_COUPLING_H__

/*============================================================================
 * Common functionnality for various coupling types.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Synchronize with applications in the same PLE coupling group.
 *
 * This function should be called before starting a new time step. The
 * current time step id is that of the last finished time step, or 0 at
 * initialization.
 *
 * Fortran Interface:
 *
 * subroutine cplsyn (ntcmabs, ntcabs, dtref)
 * *****************
 *
 * integer          ntmabs      : <-> : maximum iteration number
 * integer          ntcabs      : <-- : current iteration number
 * double precision dtref       : <-> : reference time step value
 *----------------------------------------------------------------------------*/

void CS_PROCF(cplsyn, CPLSYN)
(
 cs_int_t         *ntmabs,
 const cs_int_t   *ntcabs,
 cs_real_t        *dtref
 );

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI root communicator.
 *
 * parameters:
 *   app_name <-- name of this instance of Code_Saturne.
 *----------------------------------------------------------------------------*/

void
cs_coupling_discover_mpi_apps(const char  *app_name);

/*----------------------------------------------------------------------------
 * Finalize MPI coupling helper structures.
 *----------------------------------------------------------------------------*/

void
cs_coupling_finalize(void);

/*----------------------------------------------------------------------------
 * Return info on other applications in the same MPI root communicator.
 *
 * returns:
 *   info on other applications structure.
 *----------------------------------------------------------------------------*/

const ple_coupling_mpi_set_t *
cs_coupling_get_mpi_apps(void);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Synchronize with applications in the same PLE coupling group.
 *
 * This function should be called before starting a new time step. The
 * current time step id is that of the last finished time step, or 0 at
 * initialization.
 *
 * Default synchronization flags indicating a new iteration or end of
 * calculation are set automatically, but the user may set additional flags
 * to this function if necessary.
 *
 * parameters:
 *   flags         <-- optional additional synchronization flags
 *   current_ts_id <-- current time step id
 *   max_ts_id     <-> maximum time step id
 *   ts            <-> suggested time step value
 *----------------------------------------------------------------------------*/

void
cs_coupling_sync_apps(int      flags,
                      int      current_ts_id,
                      int     *max_ts_id,
                      double  *ts);

/*----------------------------------------------------------------------------
 * Compute extents of a mesh representation
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with mesh:
 *                     x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

ple_lnum_t
cs_coupling_mesh_extents(const void  *mesh,
                         ple_lnum_t   n_max_extents,
                         double       tolerance,
                         double       extents[]);

/*----------------------------------------------------------------------------
 * Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to the id of a given element + 1 in
 * concatenated sections of same element dimension.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   tolerance    <-- associated tolerance
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, 0 - 1 if inside,
 *                    and > 1 if outside a volume element, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh(const void         *mesh,
                          double              tolerance,
                          ple_lnum_t          n_points,
                          const ple_coord_t   point_coords[],
                          ple_lnum_t          location[],
                          float               distance[]);

/*----------------------------------------------------------------------------
 * Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to parent element numbers.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   tolerance    <-- associated tolerance
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, 0 - 1 if inside,
 *                    and > 1 if outside a volume element, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh_p(const void         *mesh,
                            double              tolerance,
                            ple_lnum_t          n_points,
                            const ple_coord_t   point_coords[],
                            ple_lnum_t          location[],
                            float               distance[]);

/*----------------------------------------------------------------------------
 * Find elements in a given mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * This function currently only handles elements of lower dimension than
 * the spatial dimension.
 *
 * Location is relative to the id of a given element + 1 in
 * concatenated sections of same element dimension.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_closest_mesh(const void         *mesh,
                               ple_lnum_t          n_points,
                               const ple_coord_t   point_coords[],
                               ple_lnum_t          location[],
                               float               distance[]);

/*----------------------------------------------------------------------------
 * Find elements in a given mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * This function currently only handles elements of lower dimension than
 * the spatial dimension.
 *
 * Location is relative to parent element numbers.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_closest_mesh_p(const void         *mesh,
                                 ple_lnum_t          n_points,
                                 const ple_coord_t   point_coords[],
                                 ple_lnum_t          location[],
                                 float               distance[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COUPLING_H__ */
