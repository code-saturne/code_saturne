#ifndef __PLE_LOCATOR_H__
#define __PLE_LOCATOR_H__

/*============================================================================
 * Locate points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2018  EDF S.A.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

#if defined(PLE_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* PLE option types */

typedef enum {

  PLE_LOCATOR_NUMBERING,
  PLE_LOCATOR_N_OPTIONS

} ple_locator_option_t;

/*----------------------------------------------------------------------------
 * Query number of extents and compute extents of a mesh representation.
 *
 * For future optimizations, computation of extents should not be limited
 * to mesh extents, but to 1 to n extents, allowing different extent
 * refinements, from global mesh to individual element extents.
 *
 * The minimum required functionality for this function is to compute
 * whole mesh extents, but it could also return extents of individual
 * elements, or intermediate extents of mesh subdivisions or coarsened
 * elements. As such, it takes an argument indicating the maximum
 * local number of extents it should compute (based on the size of
 * the extents array argument), but returns the number of extents
 * really computed, which may be lower (usually 1 for mesh extents,
 * possibly even 0 if the local mesh is empty). If n_max_extents = 1,
 * the whole mesh extents should be computed.
 *
 * If n_max_extents is set to a negative value (-1), no extents are computed,
 * but the function returns the maximum number of extents it may compute.
 * This query mode allows for the caller to allocate the correct amount
 * of memory for a subsequent call.
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with the mesh or elements (or even
 *                     aggregated elements in case of coarser representation):
 *                     x_min_0, y_min_0, ..., x_max_i, y_max_i, ...
 *                     (size: 2*dim*n_max_extents), ignored in query mode
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

typedef ple_lnum_t
(ple_mesh_extents_t) (const void  *mesh,
                      ple_lnum_t   n_max_extents,
                      double       tolerance,
                      double       extents[]);

/*----------------------------------------------------------------------------
 * Find elements in a given local mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_nodal         <-- pointer to nodal mesh representation structure
 *   tolerance_base     <-- associated base tolerance (used for bounding
 *                          box check only, not for location test)
 *   tolerance_fraction <-- associated fraction of element bounding boxes
 *                          added to tolerance
 *   n_points           <-- number of points to locate
 *   point_coords       <-- point coordinates (interleaved)
 *   point_tag          <-- optional point tag (size: n_points)
 *   location           <-> number of element containing or closest to each
 *                          point (size: n_points)
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          and > 1 if outside a volume element, or absolute
 *                          distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

typedef void
(ple_mesh_elements_locate_t) (const void         *mesh,
                              float               tolerance_base,
                              float               tolerance_fraction,
                              ple_lnum_t          n_points,
                              const ple_coord_t   point_coords[],
                              const int           point_tag[],
                              ple_lnum_t          location[],
                              float               distance[]);

/*----------------------------------------------------------------------------
 * Function pointer type for user definable logging/profiling type functions
 *----------------------------------------------------------------------------*/

typedef int
(ple_locator_log_t) (int         event,
                     int         data,
                     const char *string);

/*----------------------------------------------------------------------------
 * Structure defining a locator
 *----------------------------------------------------------------------------*/

typedef struct _ple_locator_t ple_locator_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a locator structure.
 *
 * Note that depending on the choice of ranks of the associated communicator,
 * distant ranks may in fact be truly distant or not. If n_ranks = 1 and
 * start_rank is equal to the current rank in the communicator, the locator
 * will work only locally.
 *
 * parameters:
 *   comm       <-- associated MPI communicator
 *   n_ranks    <-- number of MPI ranks associated with distant location
 *   start_rank <-- first MPI rank associated with distant location
 *
 * returns:
 *   pointer to locator
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)

ple_locator_t *
ple_locator_create(MPI_Comm  comm,
                   int       n_ranks,
                   int       start_rank);

#else

ple_locator_t *
ple_locator_create(void);

#endif

/*----------------------------------------------------------------------------
 * Destruction of a locator structure.
 *
 * parameters:
 *   this_locator <-> locator to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

ple_locator_t *
ple_locator_destroy(ple_locator_t  * this_locator);

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given mesh representation.
 *
 * parameters:
 *   this_locator       <-> pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   options            <-- options array (size PLE_LOCATOR_N_OPTIONS),
 *                          or NULL
 *   tolerance_base     <-- associated base tolerance (used for bounding
 *                          box check only, not for location test)
 *   tolerance_fraction <-- associated fraction of element bounding boxes
 *                          added to tolerance
 *   dim                <-- spatial dimension of mesh and points to locate
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_tag          <-- optional point tag (size: n_points)
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   distance           --> optional distance from point to matching element:
 *                          < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                          outside a volume element, or absolute distance
 *                          to a surface element (size: n_points)
 *   mesh_extents_f     <-- pointer to function computing mesh extents
 *   locate_f           <-- pointer to function wich updates the location[]
 *                          and distance[] arrays associated with a set of
 *                          points for points that are in an element of this
 *                          mesh, or closer to one than to previously
 *                          encountered elements.
 *----------------------------------------------------------------------------*/

void
ple_locator_set_mesh(ple_locator_t               *this_locator,
                     const void                  *mesh,
                     const int                   *options,
                     float                        tolerance_base,
                     float                        tolerance_fraction,
                     int                          dim,
                     ple_lnum_t                   n_points,
                     const ple_lnum_t             point_list[],
                     const int                    point_tag[],
                     const ple_coord_t            point_coords[],
                     float                        distance[],
                     ple_mesh_extents_t          *mesh_extents_f,
                     ple_mesh_elements_locate_t  *mesh_elements_locate_f);

/*----------------------------------------------------------------------------
 * Extend search for a locator for which set_mesh has already been called.
 *
 * parameters:
 *   this_locator       <-> pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   options            <-- options array (size PLE_LOCATOR_N_OPTIONS),
 *                          or NULL
 *   tolerance_base     <-- associated base tolerance (used for bounding
 *                          box check only, not for location test)
 *   tolerance_fraction <-- associated fraction of element bounding boxes
 *                          added to tolerance
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_tag          <-- optional point tag (size: n_points)
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   distance           --> optional distance from point to matching element:
 *                          < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                          outside a volume element, or absolute distance
 *                          to a surface element (size: n_points)
 *   mesh_extents_f     <-- pointer to function computing mesh extents
 *   locate_f           <-- pointer to function wich updates the location[]
 *                          and distance[] arrays associated with a set of
 *                          points for points that are in an element of this
 *                          mesh, or closer to one than to previously
 *                          encountered elements.
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_extend_search(ple_locator_t               *this_locator,
                          const void                  *mesh,
                          const int                   *options,
                          float                        tolerance_base,
                          float                        tolerance_fraction,
                          ple_lnum_t                   n_points,
                          const ple_lnum_t             point_list[],
                          const ple_lnum_t             point_tag[],
                          const ple_coord_t            point_coords[],
                          float                        distance[],
                          ple_mesh_extents_t          *mesh_extents_f,
                          ple_mesh_elements_locate_t  *mesh_locate_f);

/*----------------------------------------------------------------------------
 * Shift location ids for located points after locator initialization.
 *
 * This is useful mainly to switch between 0-based to 1-based numberings.
 *
 * parameters:
 *   this_locator   <-> pointer to locator structure
 *   location_shift <-- shift value
 *----------------------------------------------------------------------------*/

void
ple_locator_shift_locations(ple_locator_t  *this_locator,
                            ple_lnum_t      location_shift);

/*----------------------------------------------------------------------------
 * Return number of distant points after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of distant points.
 *----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_dist_points(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return an array of local element numbers containing (or nearest to)
 * each distant point after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   local element numbers associated with distant points.
 *----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_dist_locations(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return an array of coordinates of each distant point after
 * locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   coordinate array associated with distant points (interlaced).
 *----------------------------------------------------------------------------*/

const ple_coord_t *
ple_locator_get_dist_coords(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return number of points located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points located.
 *----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_interior(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return list of points located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points located.
 *----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_interior_list(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return number of points not located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points not located.
 *----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_exterior(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points not located.
 *----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_exterior_list(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Discard list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *----------------------------------------------------------------------------*/

void
ple_locator_discard_exterior(ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * The caller should have defined the values of distant_var[] for the
 * distant points, whose coordinates are given by
 * ple_locator_get_dist_coords(), and which are located in the elements
 * whose numbers are given by ple_locator_get_dist_locations().
 *
 * The local_var[] is defined at the located points (those whose
 * numbers are returned by ple_locator_get_interior_list().
 *
 * If the optional local_list indirection is used, it is assumed to use
 * the same base numbering as that defined by the options for the previous
 * call to ple_locator_set_mesh() or ple_locator_extend_search().
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *                     size: n_dist_points*stride
 *   local_var     <-> variable defined on located local points (received)
 *                     size: n_interior*stride
 *   local_list    <-- optional indirection list for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

void
ple_locator_exchange_point_var(ple_locator_t     *this_locator,
                               void              *distant_var,
                               void              *local_var,
                               const ple_lnum_t  *local_list,
                               size_t             type_size,
                               size_t             stride,
                               int                reverse);

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * In parallel mode, this includes communication time.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
ple_locator_get_times(const ple_locator_t  *this_locator,
                      double               *location_wtime,
                      double               *location_cpu_time,
                      double               *exchange_wtime,
                      double               *exchange_cpu_time);

/*----------------------------------------------------------------------------
 * Return communication timing information.
 *
 * In serial mode, returned times are always zero..
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
ple_locator_get_comm_times(const ple_locator_t  *this_locator,
                           double               *location_wtime,
                           double               *location_cpu_time,
                           double               *exchange_wtime,
                           double               *exchange_cpu_time);

/*----------------------------------------------------------------------------
 * Dump printout of a locator structure.
 *
 * parameters:
 *   this_locator  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
ple_locator_dump(const ple_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Get the maximum number of exchanging ranks for which we use asynchronous
 * MPI sends and receives instead of MPI_SendRecv.
 *
 * returns:
 *   the maximum number of ranks allowing asynchronous exchanges
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)

int
ple_locator_get_async_threshold(void);

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Set the maximum number of exchanging ranks for which we use asynchronous
 * MPI sends and receives instead of MPI_SendRecv.
 *
 * parameters:
 *   threshold  <-- maximum number of ranks allowing asynchronous exchanges
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)

void
ple_locator_set_async_threshold(int threshold);

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Register communication logging functions for locator instrumentation.
 *
 * By default, locators are not instrumented.
 *
 * Functions using MPE may be defined and used, but other similar systems
 * may be used.
 *
 * parameters:
 *   log_function  <-- pointer to logging function
 *   start_p_comm  <-- point to point communication start event number
 *   end_p_comm    <-- point to point communication end event number
 *   start_g_comm  <-- global communication start event number
 *   end_g_comm    <-- global communication end event number
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)

void
ple_locator_set_comm_log(ple_locator_log_t  *log_function,
                         int                 start_p_comm,
                         int                 end_p_comm,
                         int                 start_g_comm,
                         int                 end_g_comm);

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_LOCATOR_H__ */
