#ifndef __CS_COUPLING_H__
#define __CS_COUPLING_H__

/*============================================================================
 * Common functionnality for various coupling types.
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

/*----------------------------------------------------------------------------
 * Indicate if there are synchronized applications in the same
 * PLE coupling group.
 *
 * Fortran Interface:
 *
 * subroutine cplact (isync)
 * *****************
 *
 * integer          isync       : <-- : 1 if synchronized, 0 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF(cplact, CPLACT)
(
 cs_int_t         *isync
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
cs_coupling_discover_mpi_apps(const char  *app_name,
                              const char  *forced_app_type);

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
 * Return the optional synchronization flag for external couplings.
 *
 * See cs_coupling_set_sync_flag() for details.
 *
 * returns:
 *   synchronization flag to apply to couplings
 *----------------------------------------------------------------------------*/

int
cs_coupling_get_sync_flag(void);

/*----------------------------------------------------------------------------
 * Define an optional synchronization flag for external couplings.
 *
 * This flag is used by all couplings based on the PLE (Parallel Location
 * and Exchange) group synchronization mechanism, which include couplings
 * with SYRTHES 4, Code_Saturne, and NEPTUNE_CFD.
 *
 * It is defined by a mask, so for example flags f1, f2, and f3 may be
 * combined using the "f1 | f2 | f2" syntax.
 *
 * Note also that for Code_Saturne, in the case of a variable time step,
 * the reference time step is synchronized at the beginning of each
 * iteration, but the actual time step is recomputed later.
 *
 * Possible flags are:
 *   PLE_COUPLING_TS_MIN        Use smallest time step
 *   PLE_COUPLING_TS_LEADER     Prescribe time step for the group
 *                              (only one member may set this flag)
 *   PLE_COUPLING_UNSTEADY      Inform others that this instance is
 *                              using an unsteady solution approach
 *   PLE_COUPLING_STEADY        Inform others that this instance is
 *                              using a teady solution approach
 *   PLE_COUPLING_USER_1        User definable flag
 *   PLE_COUPLING_USER_2        User definable flag
 *   PLE_COUPLING_USER_3        User definable flag
 *   PLE_COUPLING_USER_4        User definable flag
 *
 * To force stopping, PLE_COUPLING_STOP may be set. In this case,
 * the calculation will stop at the first synchronization, even if
 * this function is called again with another flag.
 *
 * parameters:
 *   flag <-- synchronization flag to apply to couplings
 *----------------------------------------------------------------------------*/

void
cs_coupling_set_sync_flag(int flag);

/*----------------------------------------------------------------------------
 * Return the time step multiplier for external couplings.
 *
 * See cs_coupling_get_ts_multiplier() for details.
 *
 * returns:
 *   time step multiplier for external couplings
 *----------------------------------------------------------------------------*/

double
cs_coupling_get_ts_multiplier(void);

/*----------------------------------------------------------------------------
 * Define a time step multiplier for external couplings.
 *
 * The apparent time step for the current instance times (as viewed by
 * coupled codes) is equal to the true time step times this multiplier.
 *
 * If the synchronization flag contains "time step min" (PLE_COUPLING_TS_MIN),
 * the apparent time step is used to determine which code has the smallest
 * time step.
 *
 * parameters:
 *   m <-- time step multipier to aply to couplings
 *----------------------------------------------------------------------------*/

void
cs_coupling_set_ts_multiplier(double m);

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
 * Indicate is synchronization with applications in the same
 * PLE group is active.
 *
 * return:
 *   true if synchronization is required, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_coupling_is_sync_active(void);

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
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated base tolerance (used for bounding
 *                          box check only, not for location test)
 *   tolerance_fraction <-- associated fraction of element bounding boxes
 *                          added to tolerance
 *   n_points           <-- number of points to locate
 *   point_coords       <-- point coordinates
 *   point_tag          <-- optional point tag (size: n_points)
 *   location           <-> number of element containing or closest to each
 *                          point (size: n_points)
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          and > 1 if outside a volume element, or absolute
 *                          distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh(const void         *mesh,
                          float               tolerance_base,
                          float               tolerance_fraction,
                          ple_lnum_t          n_points,
                          const ple_coord_t   point_coords[],
                          const int           point_tag[],
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
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated base tolerance (used for bounding
 *                          box check only, not for location test)
 *   tolerance_fraction <-- associated fraction of element bounding boxes
 *                          added to tolerance
 *   n_points           <-- number of points to locate
 *   point_coords       <-- point coordinates
 *   point_tag          <-- optional point tag (size: n_points)
 *   location           <-> number of element containing or closest to each
 *                          point (size: n_points)
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          and > 1 if outside a volume element, or absolute
 *                          distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh_p(const void         *mesh,
                            float               tolerance_base,
                            float               tolerance_fraction,
                            ple_lnum_t          n_points,
                            const ple_coord_t   point_coords[],
                            const int           point_tag[],
                            ple_lnum_t          location[],
                            float               distance[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COUPLING_H__ */
