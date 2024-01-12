/*============================================================================
 * Locate points in a representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2024  EDF S.A.

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

/*!
 * \file ple_locator.c
 *
 * \brief Locate points in a representation associated with a mesh.
 */

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (__STDC_VERSION__ <202311L)
# include <stdbool.h>
#endif

#if defined(PLE_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_config_defs.h"
#include "ple_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "ple_locator.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/* Geometric operation macros*/

#define _MODULE(vect) \
  sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2])

/*
 * Algorithm ID's
 */

#define _LOCATE_BB_SENDRECV          100  /* bounding boxes + send-receive */
#define _LOCATE_BB_SENDRECV_ORDERED  200  /* bounding boxes + send-receive
                                             with communication ordering */

#define _EXCHANGE_SENDRECV       100  /* Sendrecv */
#define _EXCHANGE_ISEND_IRECV    200  /* Isend/Irecv/Waitall */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  int       n;       /* Number of intersecting distant ranks */
  int      *rank;    /* List of intersecting distant ranks */
  double   *extents; /* List of intersecting distant extents */

} _rank_intersects_t;

/*----------------------------------------------------------------------------
 * Structure defining a locator
 *----------------------------------------------------------------------------*/

struct _ple_locator_t {

  /* Basic information */
  /*-------------------*/

  int       dim;                 /* Spatial dimension */
  int       have_tags;           /* Do we use tags */

  int       locate_algorithm;    /* Location algorithm id */
  int       exchange_algorithm;  /* Exchange algorithm id */
  int       async_threshold;     /* Threshold for asynchronous exchange */

#if defined(PLE_HAVE_MPI)
  MPI_Comm  comm;                /* Associated MPI communicator */
#endif

  int       n_ranks;             /* Number of MPI ranks of distant location */
  int       start_rank;          /* First MPI rank of distant location */

  int       n_intersects;        /* Number of intersecting distant ranks */
  int      *intersect_rank;      /* List of intersecting distant ranks */
  int      *comm_order;          /* Optional communication ordering */

  ple_lnum_t    point_id_base;   /* base numbering for (external) point ids */

  ple_lnum_t   *local_points_idx;   /* Start index of local points per rank
                                       (size: n_intersects + 1)*/
  ple_lnum_t   *distant_points_idx; /* Start index of distant points per rank
                                       (size: n_intersects + 1)*/

  ple_lnum_t   *local_point_ids;        /* Local point index for data received
                                           (with blocs starting at
                                           local_points_idx[] indexes,
                                           0 to n-1 numbering) */

  ple_lnum_t   *distant_point_location; /* Location of distant points by parent
                                           element number (with blocs starting
                                           at distant_points_idx[] indexes) */
  ple_coord_t  *distant_point_coords;   /* Coordinates of distant points
                                           (with blocs starting at
                                           distant_points_idx[]*dim indexes) */

  ple_lnum_t    n_interior;         /* Number of local points located */
  ple_lnum_t   *interior_list;      /* List of points located */
  ple_lnum_t    n_exterior;         /* Number of local points not located */
  ple_lnum_t   *exterior_list;      /* List of points not located */

  /* Timing information (2 fields/time; 0: total; 1: communication) */

  double  location_wtime[2];       /* Location Wall-clock time */
  double  location_cpu_time[2];    /* Location CPU time */
  double  exchange_wtime[2];       /* Variable exchange Wall-clock time */
  double  exchange_cpu_time[2];    /* Variable exchange CPU time */
};

/*============================================================================
 * Local function pointer type documentation
 *============================================================================*/

#ifdef DOXYGEN_ONLY

/*!
 * \brief Query number of extents and compute extents of a mesh representation.
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
 * \param[in] mesh          pointer to mesh representation structure
 * \param[in] n_max_extents maximum number of sub-extents (such as element
 *                          extents) to compute, or -1 to query
 * \param[in] tolerance     addition to local extents of each element:
 *                          extent = base_extent * (1 + tolerance)
 * \param[in] extents       extents associated with the mesh or elements (or
 *                          even aggregated elements in case of coarser
 *                          representation):
 *                          x_min_0, y_min_0, ..., x_max_i, y_max_i, ...
 *                          (size: 2*dim*n_max_extents), ignored in query mode
 *
 * \return the number of extents computed
 */

typedef ple_lnum_t
(ple_mesh_extents_t) (const void  *mesh,
                      ple_lnum_t   n_max_extents,
                      double       tolerance,
                      double       extents[]);

/*!
 * \brief Find elements in a given local mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * \param[in]      this_nodal          pointer to nodal mesh representation
 *                                     structure
 * \param[in]      tolerance_base      associated fixed tolerance
 * \param[in]      tolerance_fraction  associated fraction of element bounding
 *                                     boxes added to tolerance
 * \param[in]      n_points            number of points to locate
 * \param[in]      point_coords        point coordinates
 * \param[in, out] location            number of element containing or closest
 *                                     to each point (size: n_points)
 * \param[in, out] distance            distance from point to element indicated
 *                                     by location[]: < 0 if unlocated, >= 0
 *                                     if inside; the choice of distance metric
 *                                     is left to the calling code
 *                                     (size: n_points)
 */

typedef void
(ple_mesh_elements_locate_t) (const void         *mesh,
                              float               tolerance_base,
                              float               tolerance_fraction,
                              ple_lnum_t          n_points,
                              const ple_coord_t   point_coords[],
                              ple_lnum_t          location[],
                              float               distance[]);

/*!
 * \brief Function pointer type for user definable logging/profiling
 * type functions
 */

typedef int
(ple_locator_log_t) (int         event,
                     int         data,
                     const char *string);

#endif /* DOXYGEN_ONLY */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default location algorithm */

static int _ple_locator_location_algorithm = _LOCATE_BB_SENDRECV_ORDERED;

/* maximum number of exchanging ranks for which we use asynchronous
   MPI sends and receives instead of MPI_SendRecv */

static int _ple_locator_async_threshold = 128;

#if defined(PLE_HAVE_MPI)

/* global logging function */

static ple_locator_log_t   *_ple_locator_log_func = NULL;

/* global variables associated with communication logging */

static int  _ple_locator_log_start_p_comm = 0;
static int  _ple_locator_log_end_p_comm = 0;
static int  _ple_locator_log_start_g_comm = 0;
static int  _ple_locator_log_end_g_comm = 0;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set options for a given locator.
 *
 * \param[in]  this_locator  pointer to locator structure
 * \param[in]  key           option name
 * \param[in]  value         associated value
 */
/*----------------------------------------------------------------------------*/

static void
_parse_locator_option(const char   *key,
                      const char   *value,
                      int          *location_algorithm,
                      int          *async_threshold)
{
  if (key == NULL || value == NULL)
    return;

  if (strncmp(key, "location_algorithm", 64) == 0) {
    if (strncmp(value, "bounding_box_sendreceive", 64) == 0)
      *location_algorithm = _LOCATE_BB_SENDRECV;
    else if (strncmp(value, "bounding_box_sendreceive_ordered", 64) == 0)
      *location_algorithm = _LOCATE_BB_SENDRECV_ORDERED;
  }

  else if (strncmp(key, "exchange_async_threshold", 64) == 0) {
    int i_value = 0;
    if (sscanf(value, "%d", &i_value) == 1)
      *async_threshold = i_value;
  }
}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log communication start.
 *
 * parameters:
 *   start_p_comm <-- event number for the start of the "send/recv"
 *                    or "global send/recv" state
 *   timing       <-> 0: wall-clock total; 1: CPU total;
 *                    2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_start_comm(int      start_p_comm,
                          double   timing[4])
{
  timing[2] = ple_timer_wtime();
  timing[3] = ple_timer_cpu_time();

  if(_ple_locator_log_func != NULL)
    _ple_locator_log_func(start_p_comm, 0, NULL);
}

/*----------------------------------------------------------------------------
 * Log communication end.
 *
 * parameters:
 *   end_p_comm  <-- event number for the end of the "send/recv" or
 *                   "global send/recv" state
 *   timing      <-> 0: wall-clock total; 1 CPU total;
 *                   2: wall-clock timer start; 3: CPU timer start
 *----------------------------------------------------------------------------*/

inline static void
_locator_trace_end_comm(int      end_p_comm,
                        double   timing[4])
{
  if(_ple_locator_log_func != NULL)
    _ple_locator_log_func(end_p_comm, 0, NULL);

  timing[0] += ple_timer_wtime() - timing[2];
  timing[1] += ple_timer_cpu_time() - timing[3];
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of an integer array.
 *
 * parameters:
 *   number  <-- pointer to numbers of entities that should be ordered.
 *               (if NULL, a default 1 to n numbering is considered)
 *   level   <-- level of the binary tree to descend
 *   n       <-- number of entities in the binary tree to descend
 *   order   <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_int_descend_tree(const int     number[],
                        size_t        level,
                        const size_t  n,
                        int           order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (number[i1] > number[i2]) lv_cur++;
    }

    if (lv_cur >= n) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (number[i1] >= number[i2]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of local numbers.
 *
 * parameters:
 *   number  <-- array of entity numbers (if NULL, a default 1 to n
 *               numbering is considered)
 *   order   <-- pre-allocated ordering table
 *   n       <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_int(const int     number[],
           int           order[],
           const size_t  n)
{
  size_t i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0; i < n; i++)
    order[i] = i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2);
  do {
    i--;
    _order_int_descend_tree(number, i, n, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_int_descend_tree(number, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Order communicating ranks to reduce serialization
 *
 * parameters:
 *   rank_id       <--  local rank id
 *   n_ranks       <--  communicator size
 *   n_comm_ranks  <--  number of communicating ranks
 *   comm_rank     <--  communicating ranks
 *   order         -->  communication order
 *----------------------------------------------------------------------------*/

static void
_order_comm_ranks(int        rank_id,
                  int        n_ranks,
                  int        n_comm_ranks,
                  const int  comm_rank[],
                  int        order[])
{
  int *sub_rank, *rank_dist;

  PLE_MALLOC(sub_rank, n_comm_ranks, int);
  PLE_MALLOC(rank_dist, n_comm_ranks, int);

  /* Compute ordering ("distance") key */

  for (int i = 0; i < n_comm_ranks; i++) {
    sub_rank[i] = comm_rank[i];
    rank_dist[i] = 0;
  }

  int l_r_id = rank_id;
  int m_r_id = (n_ranks+1) / 2;

  while (true) {

    if (l_r_id < m_r_id) {
      for (int i = 0; i < n_comm_ranks; i++) {
        if (sub_rank[i] >= m_r_id) {
          rank_dist[i] = rank_dist[i]*2 + 1;
          sub_rank[i] -= m_r_id;
        }
        else
          rank_dist[i] = rank_dist[i]*2;
      }
    }
    else {
      for (int i = 0; i < n_comm_ranks; i++) {
        if (sub_rank[i] >= m_r_id) {
          rank_dist[i] = rank_dist[i]*2;
          sub_rank[i] -= m_r_id;
        }
        else
          rank_dist[i] = rank_dist[i]*2 + 1;
      }
      l_r_id -= m_r_id;
    }

    if (m_r_id < 2)
      break;

    m_r_id = (m_r_id+1) / 2;

  }

  /* Order results */

  _order_int(rank_dist, order, n_comm_ranks);

  PLE_FREE(sub_rank);
  PLE_FREE(rank_dist);
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static bool
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static bool
_within_extents(int                dim,
                const ple_coord_t  coords[],
                const double       extents[])
{
  int i;
  bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim             <-- space dimension of points to locate
 *   point_list_base <-- base numbering for point list
 *   n_points        <-- number of points to locate
 *   point_list      <-- optional indirection array to point_coords
 *   point_coords    <-- coordinates of points to locate
 *                       (dimension: dim * n_points)
 *   location        <-- existing location_information, or NULL
 *   extents         --> extents associated with mesh:
 *                       x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(int                  dim,
               ple_lnum_t           point_list_base,
               ple_lnum_t           n_points,
               const ple_lnum_t     point_list[],
               const ple_coord_t    point_coords[],
               const ple_lnum_t     location[],
               double               extents[])
{
  int i;
  ple_lnum_t j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (location != NULL) {

    if (point_list != NULL) {

      for (j = 0; j < n_points; j++) {
        coord_idx = point_list[j] - point_list_base;
        for (i = 0; i < dim; i++) {
          if (extents[i]       > point_coords[(coord_idx * dim) + i])
            extents[i]       = point_coords[(coord_idx * dim) + i];
          if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
            extents[i + dim] = point_coords[(coord_idx * dim) + i];
        }
      }
    }

    else {

      for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
        for (i = 0; i < dim; i++) {
          if (extents[i]       > point_coords[(coord_idx * dim) + i])
            extents[i]       = point_coords[(coord_idx * dim) + i];
          if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
            extents[i + dim] = point_coords[(coord_idx * dim) + i];
        }
      }
    }

  }
  else {

    if (point_list != NULL) {

      for (j = 0; j < n_points; j++) {
        coord_idx = point_list[j] - point_list_base;
        for (i = 0; i < dim; i++) {
          if (extents[i]       > point_coords[(coord_idx * dim) + i])
            extents[i]       = point_coords[(coord_idx * dim) + i];
          if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
            extents[i + dim] = point_coords[(coord_idx * dim) + i];
        }
      }
    }

    else {

      for (j = 0; j < n_points; j++) {
        for (i = 0; i < dim; i++) {
          if (extents[i]       > point_coords[(j * dim) + i])
            extents[i]       = point_coords[(j * dim) + i];
          if (extents[i + dim] < point_coords[(j * dim) + i])
            extents[i + dim] = point_coords[(j * dim) + i];
        }
      }
    }

  }

}

/*----------------------------------------------------------------------------
 * Clear previous location information.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *----------------------------------------------------------------------------*/

static void
_clear_location_info(ple_locator_t  *this_locator)
{
  this_locator->n_intersects = 0;

  PLE_FREE(this_locator->intersect_rank);

  PLE_FREE(this_locator->local_points_idx);
  PLE_FREE(this_locator->distant_points_idx);

  PLE_FREE(this_locator->local_point_ids);

  PLE_FREE(this_locator->distant_point_location);
  PLE_FREE(this_locator->distant_point_coords);

  PLE_FREE(this_locator->interior_list);
  PLE_FREE(this_locator->exterior_list);
}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Update intersection rank information once location is done.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   n_points          <-- number of points to locate
 *   location_rank_id  <-> rank on which a point is located in, id
 *                         in updated locator intersection rank info out
 *----------------------------------------------------------------------------*/

static void
_location_ranks(ple_locator_t  *this_locator,
                ple_lnum_t      n_points,
                ple_lnum_t      location_rank_id[])
{
  int i, k;
  ple_lnum_t j;

  int loc_vals[2], max_vals[2];
  int comm_size;
  int *send_flag = NULL, *recv_flag = NULL, *intersect_rank_id = NULL;

  double comm_timing[4] = {0., 0., 0., 0.};

  /* Initialization */

  MPI_Comm_size(this_locator->comm, &comm_size);

  PLE_MALLOC(send_flag, comm_size, int);
  PLE_MALLOC(recv_flag, comm_size, int);

  for (i = 0; i < comm_size; i++)
    send_flag[i] = 0;

  for (j = 0; j < n_points; j++) {
    if (location_rank_id[j] > -1)
      send_flag[location_rank_id[j]] = 1;
  }

  /* As exchange detection is asymetric, synchronize it */

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  MPI_Alltoall(send_flag, 1, MPI_INT, recv_flag, 1, MPI_INT,
               this_locator->comm);

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  /* Update number of "intersects" and matching rank info */

  this_locator->n_intersects = 0;

  for (i = 0; i < this_locator->n_ranks; i++) {
    j = i + this_locator->start_rank;
    if (send_flag[j] == 1 || recv_flag[j] == 1)
      this_locator->n_intersects++;
  }

  PLE_REALLOC(this_locator->intersect_rank,
              this_locator->n_intersects,
              int);

  for (i = 0, k = 0; i < this_locator->n_ranks; i++) {
    j = i + this_locator->start_rank;
    if (send_flag[j] == 1 || recv_flag[j] == 1)
      this_locator->intersect_rank[k++] = j;
  }

  if (this_locator->locate_algorithm == _LOCATE_BB_SENDRECV_ORDERED) {
    int rank_id;
    MPI_Comm_rank(this_locator->comm, &rank_id);
    PLE_REALLOC(this_locator->comm_order, this_locator->n_intersects, int);
    _order_comm_ranks(rank_id,
                      comm_size,
                      this_locator->n_intersects,
                      this_locator->intersect_rank,
                      this_locator->comm_order);
  }

  PLE_FREE(send_flag);
  PLE_FREE(recv_flag);

  /* Now convert location rank id to intersect rank */

  PLE_MALLOC(intersect_rank_id, this_locator->n_ranks, int);

  for (i = 0; i < this_locator->n_ranks; i++)
    intersect_rank_id[i] = -1;

  for (i = 0; i < this_locator->n_intersects; i++) {
    intersect_rank_id[  this_locator->intersect_rank[i]
                      - this_locator->start_rank] = i;
  }

  for (j = 0; j < n_points; j++) {
    k = location_rank_id[j] - this_locator->start_rank;
    if (k > -1)
      location_rank_id[j] = intersect_rank_id[k];
  }

  PLE_FREE(intersect_rank_id);

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  loc_vals[0] = this_locator->n_intersects;
  loc_vals[1] = this_locator->async_threshold;

  MPI_Allreduce(loc_vals, max_vals, 2, MPI_INT, MPI_MAX,
                this_locator->comm);

  if (max_vals[0] <= max_vals[1])
    this_locator->exchange_algorithm = _EXCHANGE_ISEND_IRECV;
  else
    this_locator->exchange_algorithm = _EXCHANGE_SENDRECV;

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Determine or update possibly intersecting ranks for unlocated elements,
 * in parallel.
 *
 * parameters:
 *   this_locator       <-- pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated fixed tolerance
 *   tolerance_fraction <-- associated fraction of element bounding
 *                          boxes added to tolerance
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   location           <-> number of distant element containing or closest
 *                          to each point, or -1 (size: n_points)
 *   mesh_extents_f     <-- pointer to function computing mesh or mesh
 *                          subset or element extents
 *
 * returns:
 *   local rank intersection info
 *----------------------------------------------------------------------------*/

static _rank_intersects_t
_intersects_distant(ple_locator_t       *this_locator,
                    const void          *mesh,
                    float                tolerance_base,
                    float                tolerance_fraction,
                    ple_lnum_t           n_points,
                    const ple_lnum_t     point_list[],
                    const ple_coord_t    point_coords[],
                    const ple_lnum_t     location[],
                    ple_mesh_extents_t  *mesh_extents_f)
{
  int stride2;
  double extents[12];

  int j;
  int stride4;
  int comm_rank, comm_size;
  int n_intersects;
  int  *intersect_rank;
  double *recvbuf;

  double comm_timing[4] = {0., 0., 0., 0.};
  const int dim = this_locator->dim;

  _rank_intersects_t intersects;

  /* Update intersects */

  intersects.n = 0;

  /* initialize mesh extents in case mesh is empty or dim < 3 */
  for (int i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  mesh_extents_f(mesh,
                 1,
                 tolerance_fraction,
                 extents);

  _point_extents(dim,
                 this_locator->point_id_base,
                 n_points,
                 point_list,
                 point_coords,
                 location,
                 extents + 2*dim);

  for (int i = 0; i < dim; i++) {

    if (extents[i] > -HUGE_VAL + tolerance_base)
      extents[i]         -= tolerance_base;
    else
      extents[i] = -HUGE_VAL;
    if (extents[i] < HUGE_VAL - tolerance_base)
      extents[i +   dim] += tolerance_base;
    else
      extents[i +   dim] = HUGE_VAL;

    if (extents[i + 2*dim] > -HUGE_VAL + tolerance_base)
      extents[i + 2*dim] -= tolerance_base;
    else
      extents[i + 2*dim] = -HUGE_VAL;
    if (extents[i + 3*dim] < HUGE_VAL - tolerance_base)
      extents[i + 3*dim] += tolerance_base;
    else
      extents[i + 3*dim] = HUGE_VAL;

  }

  /* Exchange extent information */

  MPI_Comm_rank(this_locator->comm, &comm_rank);
  MPI_Comm_size(this_locator->comm, &comm_size);

  stride2 = dim * 2; /* Stride for one type of extent */
  stride4 = dim * 4; /* Stride for element and vertex
                        extents, end-to-end */

  PLE_MALLOC(recvbuf, stride4*comm_size, double);

  _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

  MPI_Allgather(extents, stride4, MPI_DOUBLE, recvbuf, stride4, MPI_DOUBLE,
                this_locator->comm);

  _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

  /* Count and mark possible overlaps */

  n_intersects = 0;
  PLE_MALLOC(intersect_rank, this_locator->n_ranks, int);

  for (int i = 0; i < this_locator->n_ranks; i++) {
    j = this_locator->start_rank + i;
    if (  (_intersect_extents(dim,
                              extents + (2*dim),
                              recvbuf + (j*stride4)) == true)
        || (_intersect_extents(dim,
                               extents,
                               recvbuf + (j*stride4) + (2*dim)) == true)) {
      intersect_rank[n_intersects] = j;
      n_intersects += 1;
    }
  }

  intersects.n = n_intersects;
  PLE_MALLOC(intersects.rank, intersects.n, int);
  PLE_MALLOC(intersects.extents, intersects.n * stride2, double);

  for (int i = 0; i < intersects.n; i++) {

    intersects.rank[i] = intersect_rank[i];

    /* Copy only distant element (and not point) extents */

    for (j = 0; j < stride2; j++)
      intersects.extents[i*stride2 + j]
        = recvbuf[intersect_rank[i]*stride4 + j];

  }

  /* Free temporary memory */

  PLE_FREE(intersect_rank);
  PLE_FREE(recvbuf);

  /* Finalize timing */

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];

  return intersects;
}

/*----------------------------------------------------------------------------
 * Initialize location information from previous locator info,
 * in parallel mode.
 *
 * parameters:
 *   this_locator     <-> pointer to locator structure
 *   n_points         <-- number of points to locate
 *   location         --> number of distant element containing or closest
 *                        to each point, or -1 (size: n_points)
 *   location_rank_id --> rank id for distant element containing or closest
 *                        to each point, or -1
 *----------------------------------------------------------------------------*/

static void
_transfer_location_distant(ple_locator_t  *this_locator,
                           ple_lnum_t      n_points,
                           ple_lnum_t      location[],
                           ple_lnum_t      location_rank_id[])
{
  int dist_rank;
  ple_lnum_t j, k, n_points_loc, n_points_dist, dist_v_idx;

  double comm_timing[4] = {0., 0., 0., 0.};

  const ple_lnum_t idb = this_locator->point_id_base;
  const ple_lnum_t *_interior_list = this_locator->interior_list;

  /* Initialize locations */

  for (j = 0; j < n_points; j++) {
    location[j] = -1;
    location_rank_id[j] = -1;
  }

  /* Update with existing location info *
     exact location info is not necessary, so not determined;
     all that is necessary is to mark located points as such */

  for (int li = 0; li < this_locator->n_intersects; li++) {

    int i = (this_locator->comm_order != NULL) ?
      this_locator->comm_order[li] : li;

    MPI_Status status;
    ple_lnum_t *loc_v_buf, *dist_v_ptr;
    const ple_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    dist_v_idx = this_locator->distant_points_idx[i];

    PLE_MALLOC(loc_v_buf, n_points_loc, ple_lnum_t);

    /* Exchange information */

    dist_v_ptr = this_locator->distant_point_location + dist_v_idx;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(dist_v_ptr, n_points_dist, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 loc_v_buf, n_points_loc, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    for (k = 0; k < n_points_loc; k++) {
      ple_lnum_t pt_id = _interior_list[_local_point_ids[k]] - idb;
      location[pt_id] = loc_v_buf[k];
      location_rank_id[pt_id] = dist_rank;
    }

    PLE_FREE(loc_v_buf);
  } /* End of loop on MPI ranks */

  this_locator->n_intersects = 0;
  PLE_FREE(this_locator->intersect_rank);
  PLE_FREE(this_locator->comm_order);
  PLE_FREE(this_locator->local_points_idx);
  PLE_FREE(this_locator->distant_points_idx);
  PLE_FREE(this_locator->local_point_ids);
  PLE_FREE(this_locator->distant_point_location);
  PLE_FREE(this_locator->distant_point_coords);

  this_locator->n_interior = 0;
  this_locator->n_exterior = 0;
  PLE_FREE(this_locator->interior_list);
  PLE_FREE(this_locator->exterior_list);
}

/*----------------------------------------------------------------------------
 * Location of points not yet located on the closest elements.
 *
 * parameters:
 *   this_locator       <-> pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated fixed tolerance
 *   tolerance_fraction <-- associated fraction of element bounding
 *                          boxes added to tolerance
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_tag          <-- optional point tag (size: n_points)
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   location           <-> number of distant element containing or closest
 *                          to each point, or -1 (size: n_points)
 *   location_rank_id   <-> rank id for distant element containing or closest
 *                          to each point, or -1
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *   mesh_extents_f     <-- function computing mesh or mesh subset extents
 *   mesh_locate_f      <-- function locating the points on local elements
 *----------------------------------------------------------------------------*/

static void
_locate_distant(ple_locator_t               *this_locator,
                const void                  *mesh,
                float                        tolerance_base,
                float                        tolerance_fraction,
                ple_lnum_t                   n_points,
                const ple_lnum_t             point_list[],
                const int                    point_tag[],
                const ple_coord_t            point_coords[],
                ple_lnum_t                   location[],
                ple_lnum_t                   location_rank_id[],
                float                        distance[],
                ple_mesh_extents_t          *mesh_extents_f,
                ple_mesh_elements_locate_t  *mesh_locate_f)
{
  int k;
  int dist_rank, dist_index;
  ple_lnum_t j;
  ple_lnum_t n_coords_loc, n_coords_dist;
  ple_lnum_t *location_loc, *location_dist;
  int *tag_dist, *send_tag;
  ple_coord_t *coords_dist, *send_coords;
  float *distance_dist, *distance_loc;

  _rank_intersects_t intersects;
  MPI_Status status;

  ple_lnum_t _n_points = 0;
  double extents[6];

  double comm_timing[4] = {0., 0., 0., 0.};

  ple_lnum_t *_point_list = NULL, *_point_id = NULL, *send_id = NULL;
  const ple_lnum_t *_point_list_p = NULL;

  const int dim = this_locator->dim;
  const int stride = dim * 2;
  const int have_tags = this_locator->have_tags;
  const ple_lnum_t idb = this_locator->point_id_base;

  /* Filter non-located points */

  _n_points = 0;
  _point_list_p = point_list;

  for (j = 0; j < n_points; j++) {
    if (location[j] < 0)
      _n_points++;
  }

  if (_n_points < n_points) {

    PLE_MALLOC(_point_list, _n_points, ple_lnum_t);
    _point_list_p = _point_list;

    _n_points = 0;
    if (point_list == NULL) {
      _point_id = _point_list;
      for (j = 0; j < n_points; j++) {
        if (location[j] < 0)
          _point_list[_n_points++] = j + idb;
      }
    }
    else {
      PLE_MALLOC(_point_id, _n_points, ple_lnum_t);
      for (j = 0; j < n_points; j++) {
        if (location[j] < 0) {
          _point_list[_n_points] = point_list[j];
          _point_id[_n_points] = j + idb;
          _n_points++;
        }
      }
    }

  }

  /* Update intersect for current point list */

  intersects = _intersects_distant(this_locator,
                                   mesh,
                                   tolerance_base,
                                   tolerance_fraction,
                                   _n_points,
                                   _point_list_p,
                                   point_coords,
                                   location,
                                   mesh_extents_f);

  /* Allocate buffers */

  PLE_MALLOC(send_coords, _n_points * dim, ple_coord_t);
  PLE_MALLOC(send_id, _n_points, ple_lnum_t);

  if (have_tags)
    PLE_MALLOC(send_tag, _n_points, int);
  else
    send_tag = NULL;

  int *comm_order = NULL;
  if (this_locator->locate_algorithm == _LOCATE_BB_SENDRECV_ORDERED) {
    int comm_size, rank_id;
    MPI_Comm_size(this_locator->comm, &comm_size);
    MPI_Comm_rank(this_locator->comm, &rank_id);
    PLE_MALLOC(comm_order, intersects.n, int);
    _order_comm_ranks(rank_id,
                      comm_size,
                      intersects.n,
                      intersects.rank,
                      comm_order);
  }

  /* First loop on possibly intersecting distant ranks */
  /*---------------------------------------------------*/

  for (int li = 0; li < intersects.n; li++) {

    int i = (comm_order != NULL) ?
      comm_order[li] : li;

    dist_index = i; /* Ordering (communication schema) not yet optimized */
    dist_rank  = intersects.rank[dist_index];

    /* Prepare and send coords that should fit in each send buffer */
    /* Reset buffers for current intersect rank */

    n_coords_loc = 0;

    for (k = 0; k < stride; k++)
      extents[k] = intersects.extents[dist_index*stride + k];

    /* Build partial buffer */

    for (j = 0; j < _n_points; j++) {

      ple_lnum_t coord_idx;

      if (_point_list_p != NULL)
        coord_idx = _point_list_p[j] - idb;
      else
        coord_idx = j;

      if (_within_extents(dim,
                          &(point_coords[dim*coord_idx]),
                          extents) == true) {

        if (_point_id != NULL)
          send_id[n_coords_loc] = _point_id[j] -idb;
        else
          send_id[n_coords_loc] = j;

        for (k = 0; k < dim; k++)
          send_coords[n_coords_loc*dim + k] = point_coords[dim*coord_idx + k];

        if (have_tags)
          send_tag[n_coords_loc] = point_tag[j];

        n_coords_loc += 1;
      }

    }

    /* Send then receive partial buffer */

    dist_rank = intersects.rank[dist_index];

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 &n_coords_dist, 1, PLE_MPI_LNUM, dist_rank,
                 PLE_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_MALLOC(coords_dist, n_coords_dist*dim, ple_coord_t);
    if (have_tags)
      PLE_MALLOC(tag_dist, n_coords_dist, int);
    else
      tag_dist = NULL;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 coords_dist, (int)(n_coords_dist*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    if (have_tags)
      MPI_Sendrecv(send_tag, (int)(n_coords_loc),
                   MPI_INT, dist_rank, PLE_MPI_TAG,
                   tag_dist, (int)(n_coords_dist),
                   MPI_INT, dist_rank, PLE_MPI_TAG,
                   this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    /* Now locate received coords on local rank */

    PLE_MALLOC(location_dist, n_coords_dist, ple_lnum_t);
    PLE_MALLOC(distance_dist, n_coords_dist, float);

    for (j = 0; j < n_coords_dist; j++) {
      location_dist[j] = -1;
      distance_dist[j] = -1.0;
    }

    mesh_locate_f(mesh,
                  tolerance_base,
                  tolerance_fraction,
                  n_coords_dist,
                  coords_dist,
                  tag_dist,
                  location_dist,
                  distance_dist);

    PLE_FREE(tag_dist);
    PLE_FREE(coords_dist);

    /* Exchange location return information with distant rank */

    PLE_MALLOC(location_loc, n_coords_loc, ple_lnum_t);
    PLE_MALLOC(distance_loc, n_coords_loc, float);

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(location_dist, (int)n_coords_dist,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 location_loc, (int)n_coords_loc,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(distance_dist, (int)n_coords_dist,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 distance_loc, (int)n_coords_loc,
                 MPI_FLOAT, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_FREE(location_dist);
    PLE_FREE(distance_dist);

    /* Now update location information */

    for (j = 0; j < n_coords_loc; j++) {

      ple_lnum_t l = send_id[j];

      if (   (distance_loc[j] > -0.1)
          && (distance_loc[j] < distance[l] || location[l] < 0)) {
        location_rank_id[l] = dist_rank;
        location[l] = location_loc[j];
        distance[l] = distance_loc[j];
      }

    }

    PLE_FREE(location_loc);
    PLE_FREE(distance_loc);

  }

  PLE_FREE(comm_order);

  /* Free temporary arrays */

  PLE_FREE(send_id);
  PLE_FREE(send_tag);
  PLE_FREE(send_coords);

  if (_point_list != point_list) {
    if (_point_id != _point_list)
      PLE_FREE(_point_id);
    PLE_FREE(_point_list);
  }

  PLE_FREE(intersects.rank);
  PLE_FREE(intersects.extents);

  /* Finalize timing */

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given mesh representation and point set.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   mesh              <-- pointer to mesh representation structure
 *   dim               <-- spatial dimension
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *   point_tag         <-- optional point tag (size: n_points)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *   location          <-> number of distant element containing or closest
 *                          to each point, or -1 (size: n_points)
 *   location_rank_id  <-> rank id for distant element containing or closest
 *                          to each point, or -1
 *   distance          --> optional distance from point to matching element:
 *                         < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                         outside a volume element, or absolute distance
 *                         to a surface element (size: n_points)
 *   mesh_extents_f    <-- pointer to function computing mesh or mesh
 *                         subset or element extents
 *   mesh_locate_f     <-- function locating points in or on elements
 *----------------------------------------------------------------------------*/

static void
_locate_all_distant(ple_locator_t               *this_locator,
                    const void                  *mesh,
                    float                        tolerance_base,
                    float                        tolerance_fraction,
                    ple_lnum_t                   n_points,
                    const ple_lnum_t             point_list[],
                    const int                    point_tag[],
                    const ple_coord_t            point_coords[],
                    ple_lnum_t                   location[],
                    ple_lnum_t                   location_rank_id[],
                    float                        distance[],
                    ple_mesh_extents_t          *mesh_extents_f,
                    ple_mesh_elements_locate_t  *mesh_locate_f)
{
  int k;
  int dist_rank;
  ple_lnum_t j;
  ple_lnum_t n_coords_loc, n_coords_dist, n_interior, n_exterior;
  ple_lnum_t coord_idx, start_idx;
  ple_lnum_t *send_id, *send_location;
  ple_lnum_t *location_count, *location_shift;
  float *_distance = distance;

  MPI_Status status;

  double comm_timing[4] = {0., 0., 0., 0.};
  const int dim = this_locator->dim;
  const ple_lnum_t idb = this_locator->point_id_base;

  if (distance == NULL) {
    PLE_MALLOC(_distance, n_points, float);
    for (j = 0; j < n_points; j++)
      _distance[j] = -1.0;
  }

  /* Now update locations */

  _locate_distant(this_locator,
                  mesh,
                  tolerance_base,
                  tolerance_fraction,
                  n_points,
                  point_list,
                  point_tag,
                  point_coords,
                  location,
                  location_rank_id,
                  _distance,
                  mesh_extents_f,
                  mesh_locate_f);

  /* Update info on communicating ranks and matching ids */
  /*----------------------------------------------------*/

  _location_ranks(this_locator,
                  n_points,
                  location_rank_id);

  /* Reorganize location information */
  /*---------------------------------*/

  /* Now that location is done, the location[] array contains
     either -1 if a point was not located, or a local index
     (associated with the corresponding rank); the distance[] array
     is not needed anymore now that all comparisons have been done */

  if (_distance != distance)
    PLE_FREE(_distance);

  PLE_MALLOC(location_shift, this_locator->n_intersects, ple_lnum_t);
  PLE_MALLOC(location_count, this_locator->n_intersects, ple_lnum_t);

  for (int i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    if (location_rank_id[j] > -1)
      location_count[location_rank_id[j]] += 1;
    else
      n_exterior += 1;
  }

  this_locator->n_interior = n_points - n_exterior;
  PLE_MALLOC(this_locator->interior_list, this_locator->n_interior, ple_lnum_t);

  this_locator->n_exterior = n_exterior;
  PLE_MALLOC(this_locator->exterior_list, this_locator->n_exterior, ple_lnum_t);

  if (this_locator->n_intersects > 0)
    location_shift[0] = 0;
  for (int i = 1; i < this_locator->n_intersects; i++)
    location_shift[i] = location_shift[i-1] + location_count[i-1];

  for (int i = 0; i < this_locator->n_intersects; i++)
    location_count[i] = 0;

  PLE_MALLOC(send_id, n_points, ple_lnum_t);
  PLE_MALLOC(send_location, n_points, ple_lnum_t);

  /* send_id[] will now contain information for all blocks */
  for (j = 0; j < n_points; j++)
    send_id[j] = -1;

  n_interior = 0;
  n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    const int l_rank = location_rank_id[j];
    if (l_rank > -1) {
      send_id[location_shift[l_rank] + location_count[l_rank]] = j;
      location_count[l_rank] += 1;
      this_locator->interior_list[n_interior] = j + idb;
      n_interior += 1;
    }
    else {
      this_locator->exterior_list[n_exterior] = j + idb;
      n_exterior += 1;
    }
  }

  /* Second loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  /* Count and organize total number of local and distant points */

  PLE_REALLOC(this_locator->local_points_idx,
              this_locator->n_intersects + 1,
              ple_lnum_t);

  PLE_REALLOC(this_locator->distant_points_idx,
              this_locator->n_intersects + 1,
              ple_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->distant_points_idx[0] = 0;

  for (int li = 0; li < this_locator->n_intersects; li++) {

    int i = (this_locator->comm_order != NULL) ?
      this_locator->comm_order[li] : li;

    dist_rank = this_locator->intersect_rank[i];

    n_coords_loc = location_count[i];

    this_locator->local_points_idx[i+1] = n_coords_loc;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&n_coords_loc, 1, PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 &n_coords_dist, 1, PLE_MPI_LNUM, dist_rank,
                 PLE_MPI_TAG, this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    this_locator->distant_points_idx[i+1] = n_coords_dist;

  }

  /* Counts to index */

  this_locator->local_points_idx[0] = 0;
  this_locator->distant_points_idx[0] = 0;

  for (int i = 0; i < this_locator->n_intersects; i++) {
    this_locator->local_points_idx[i+1] += this_locator->local_points_idx[i];
    this_locator->distant_points_idx[i+1] += this_locator->distant_points_idx[i];
  }

  /* Third loop on possibly intersecting distant ranks */
  /*----------------------------------------------------*/

  PLE_REALLOC(this_locator->local_point_ids,
              this_locator->local_points_idx[this_locator->n_intersects],
              ple_lnum_t);

  PLE_REALLOC(this_locator->distant_point_location,
              this_locator->distant_points_idx[this_locator->n_intersects],
              ple_lnum_t);

  PLE_REALLOC(this_locator->distant_point_coords,
              this_locator->distant_points_idx[this_locator->n_intersects] * dim,
              ple_coord_t);

  for (int li = 0; li < this_locator->n_intersects; li++) {

    ple_coord_t *send_coords;

    int i = (this_locator->comm_order != NULL) ?
      this_locator->comm_order[li] : li;
    dist_rank = this_locator->intersect_rank[i];

    n_coords_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_coords_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    start_idx = this_locator->local_points_idx[i];

    PLE_MALLOC(send_coords, n_coords_loc * dim, ple_coord_t);

    for (j = 0; j < n_coords_loc; j++) {

      coord_idx = send_id[location_shift[i] + j];
      assert(coord_idx > -1);
      this_locator->local_point_ids[start_idx + j] = coord_idx;
      send_location[j] = location[coord_idx];
      if (point_list != NULL) {
        for (k = 0; k < dim; k++)
          send_coords[j*dim + k]
            = point_coords[dim*(point_list[coord_idx] - idb) + k];
      }
      else {
        for (k = 0; k < dim; k++)
          send_coords[j*dim + k]
            = point_coords[dim*coord_idx + k];
      }
    }

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(send_location, (int)n_coords_loc,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 (this_locator->distant_point_location
                  + this_locator->distant_points_idx[i]), (int)n_coords_dist,
                 PLE_MPI_LNUM, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    MPI_Sendrecv(send_coords, (int)(n_coords_loc*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 (this_locator->distant_point_coords
                  + (this_locator->distant_points_idx[i]*dim)),
                 (int)(n_coords_dist*dim),
                 PLE_MPI_COORD, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    PLE_FREE(send_coords);

  }

  PLE_FREE(location_count);
  PLE_FREE(location_shift);

  PLE_FREE(send_location);
  PLE_FREE(send_id);

  PLE_FREE(location);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Initialize location information from previous locator info,
 * in serial mode.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   n_points          <-- number of points to locate
 *   location          --> number of distant element containing or closest
 *                         to each point, or -1 (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_transfer_location_local(ple_locator_t  *this_locator,
                         ple_lnum_t      n_points,
                         ple_lnum_t      location[])
{
  ple_lnum_t j;

  /* Initialize locations */

  for (j = 0; j < n_points; j++)
    location[j] = -1;

  /* Update with existing location info *
     exact location info is not necessary, so not determined;
     all that is necessary is to mark located points as such */

  if (this_locator->n_intersects == 1) {

    const ple_lnum_t _n_points =    this_locator->local_points_idx[1]
                                  - this_locator->local_points_idx[0];

    const ple_lnum_t idb = this_locator->point_id_base;
    const ple_lnum_t *_interior_list = this_locator->interior_list;

    const ple_lnum_t *_local_point_ids = this_locator->local_point_ids;
    const ple_lnum_t *dist_v_ptr = this_locator->distant_point_location;

    for (j = 0; j < _n_points; j++) {
      ple_lnum_t k = _interior_list[_local_point_ids[j] - idb];
      location[k] = dist_v_ptr[j];
    }

  }

  this_locator->n_intersects = 0;
  PLE_FREE(this_locator->intersect_rank);
  PLE_FREE(this_locator->comm_order);
  PLE_FREE(this_locator->local_points_idx);
  PLE_FREE(this_locator->distant_points_idx);
  PLE_FREE(this_locator->local_point_ids);
  PLE_FREE(this_locator->distant_point_location);
  PLE_FREE(this_locator->distant_point_coords);

  this_locator->n_interior = 0;
  this_locator->n_exterior = 0;
  PLE_FREE(this_locator->interior_list);
  PLE_FREE(this_locator->exterior_list);
}

/*----------------------------------------------------------------------------
 * Determine or update possibly intersecting ranks for unlocated elements,
 * in parallel.
 *
 * parameters:
 *   this_locator       <-- pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated fixed tolerance
 *   tolerance_fraction <-- associated fraction of element bounding
 *                          boxes added to tolerance
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   location           <-> number of distant element containing or closest
 *                          to each point, or -1 (size: n_points)
 *   mesh_extents_f     <-- pointer to function computing mesh or mesh
 *                          subset or element extents
 *
 * returns:
 *   local rank intersection info
 *----------------------------------------------------------------------------*/

static _rank_intersects_t
_intersects_local(ple_locator_t       *this_locator,
                  const void          *mesh,
                  float                tolerance_base,
                  float                tolerance_fraction,
                  ple_lnum_t           n_points,
                  const ple_lnum_t     point_list[],
                  const ple_coord_t    point_coords[],
                  const ple_lnum_t     location[],
                  ple_mesh_extents_t  *mesh_extents_f)
{
  int i;
  int stride2;
  double extents[12];

  int j;
  int n_intersects;

  const int dim = this_locator->dim;

  _rank_intersects_t intersects;

  /* Update intersects */

  intersects.n = 0;

  mesh_extents_f(mesh,
                 1,
                 tolerance_fraction,
                 extents);

  _point_extents(dim,
                 this_locator->point_id_base,
                 n_points,
                 point_list,
                 point_coords,
                 location,
                 extents + 2*dim);

  for (i = 0; i < dim; i++) {

    if (extents[i] > -HUGE_VAL + tolerance_base)
      extents[i]         -= tolerance_base;
    else
      extents[i] = -HUGE_VAL;
    if (extents[i] < HUGE_VAL - tolerance_base)
      extents[i +   dim] += tolerance_base;
    else
      extents[i +   dim] = HUGE_VAL;

    if (extents[i + 2*dim] > -HUGE_VAL + tolerance_base)
      extents[i + 2*dim] -= tolerance_base;
    else
      extents[i + 2*dim] = -HUGE_VAL;
    if (extents[i + 3*dim] < HUGE_VAL - tolerance_base)
      extents[i + 3*dim] += tolerance_base;
    else
      extents[i + 3*dim] = HUGE_VAL;

  }

  /* Determine possible overlap */

  stride2 = dim * 2; /* Stride for one type of extent */

  n_intersects = 0;

  if (_intersect_extents(dim, extents, extents + (2*dim)) == true)
    n_intersects += 1;

  intersects.n = n_intersects;
  PLE_MALLOC(intersects.extents, intersects.n * stride2, double);

  /* Copy only element (and not point) extents */

  for (j = 0; j < stride2; j++)
    intersects.extents[j] = extents[j];

  return intersects;
}

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given mesh representation and point set.
 *
 * parameters:
 *   this_locator       <-> pointer to locator structure
 *   mesh               <-- pointer to mesh representation structure
 *   tolerance_base     <-- associated fixed tolerance
 *   tolerance_fraction <-- associated fraction of element bounding
 *                          boxes added to tolerance
 *   n_points           <-- number of points to locate
 *   point_list         <-- optional indirection array to point_coords
 *   point_tag          <-- optional point tag (size: n_points)
 *   point_coords       <-- coordinates of points to locate
 *                          (dimension: dim * n_points)
 *   distance           --> optional distance from point to matching element:
 *                          < 0 if unlocated; 0 - 1 if inside and > 1 if
 *                          outside a volume element, or absolute distance
 *                          to a surface element (size: n_points)
 *   mesh_extents_f     <-- function computing mesh or mesh subset extents
 *   mesh_locate_f      <-- function locating the closest local elements
 *----------------------------------------------------------------------------*/

static void
_locate_all_local(ple_locator_t               *this_locator,
                  const void                  *mesh,
                  float                        tolerance_base,
                  float                        tolerance_fraction,
                  ple_lnum_t                   n_points,
                  const ple_lnum_t             point_list[],
                  const int                    point_tag[],
                  const ple_coord_t            point_coords[],
                  ple_lnum_t                   location[],
                  float                        distance[],
                  ple_mesh_extents_t          *mesh_extents_f,
                  ple_mesh_elements_locate_t  *mesh_locate_f)
{
  int l;
  ple_lnum_t j, k;
  ple_lnum_t n_coords, n_interior, n_exterior, coord_idx;
  _rank_intersects_t intersects;

  const int dim = this_locator->dim;
  const int have_tags = this_locator->have_tags;
  const ple_lnum_t idb = this_locator->point_id_base;

  /* Update intersect for current point list */

  intersects = _intersects_local(this_locator,
                                 mesh,
                                 tolerance_base,
                                 tolerance_fraction,
                                 n_points,
                                 point_list,
                                 point_coords,
                                 location,
                                 mesh_extents_f);

  n_coords = 0;

  this_locator->n_intersects = intersects.n;

  /* Build partial buffer */

  if (intersects.n > 0) {

    ple_lnum_t *_location = location;
    float *_distance = distance;

    ple_lnum_t *id;
    int *tag;
    ple_coord_t *coords;

    PLE_MALLOC(coords, n_points * dim, ple_coord_t);
    PLE_MALLOC(id, n_points, ple_lnum_t);

    if (have_tags)
      PLE_MALLOC(tag, n_points, int);
    else
      tag = NULL;

    for (j = 0; j < n_points; j++) {

      if (location[j] > -1) {
        if (distance == NULL)
          continue;
        else if (distance[j] < 1)
          continue;
      }

      if (point_list != NULL)
        coord_idx = point_list[j] - idb;
      else
        coord_idx = j;

      if (_within_extents(dim,
                          &(point_coords[dim*coord_idx]),
                          intersects.extents) == true) {

        for (k = 0; k < dim; k++)
          coords[n_coords*dim + k]
            = point_coords[dim*coord_idx + k];

        id[n_coords] = j;

        if (have_tags)
          tag[n_coords] = point_tag[j];

        n_coords += 1;
      }

    }

    PLE_REALLOC(coords, n_coords * dim, ple_coord_t);
    PLE_REALLOC(id, n_coords, ple_lnum_t);
    if (have_tags)
      PLE_REALLOC(tag, n_coords, int);

    if (n_coords < n_points) {
      PLE_MALLOC(_location, n_coords, ple_lnum_t);
      for (j = 0; j < n_coords; j++)
        _location[j] = location[id[j]];
    }
    if (distance == NULL || n_coords < n_points) {
      PLE_MALLOC(_distance, n_coords, float);
      if (distance == NULL) {
        for (j = 0; j < n_coords; j++)
          _distance[j] = -1.0;
      }
      else {
        for (j = 0; j < n_coords; j++)
          _distance[j] = distance[id[j]];
      }
    }

    mesh_locate_f(mesh,
                  tolerance_base,
                  tolerance_fraction,
                  n_coords,
                  coords,
                  tag,
                  _location,
                  _distance);

    PLE_FREE(coords);

    if (n_coords < n_points) {
      for (j = 0; j < n_coords; j++) {
        if (_location[j] > -1) {
          k = id[j];
          if (distance != NULL) {
            if (distance[k] <= _distance[j])
              continue;
            else
              distance[k] = _distance[j];
          }
          location[k] = _location[j];
        }
      }
    }

    if (_location != location)
      PLE_FREE(_location);

    if (_distance != distance)
      PLE_FREE(_distance);

    PLE_FREE(tag);
    PLE_FREE(id);

  }

  PLE_FREE(intersects.extents);

  /* Reorganize location information */
  /*---------------------------------*/

  /* Now that location is done, the location[] array contains
     either -1 if a point was not located, or a local index;
     the distance[] array is not needed anymore now that all comparisons have
     been done */

  this_locator->n_interior = 0;
  this_locator->n_exterior = 0;
  for (j = 0; j < n_points; j++) {
    if (location[j] > -1)
      this_locator->n_interior += 1;
    else
      this_locator->n_exterior += 1;
  }

  PLE_MALLOC(this_locator->interior_list, this_locator->n_interior, ple_lnum_t);
  PLE_MALLOC(this_locator->exterior_list, this_locator->n_exterior, ple_lnum_t);

  /* Organize total number of "local" and "distant" points */

  PLE_REALLOC(this_locator->local_points_idx, 2, ple_lnum_t);
  PLE_REALLOC(this_locator->distant_points_idx, 2, ple_lnum_t);

  this_locator->local_points_idx[0] = 0;
  this_locator->local_points_idx[1] = this_locator->n_interior;

  this_locator->distant_points_idx[0] = 0;
  this_locator->distant_points_idx[1] = this_locator->n_interior;

  this_locator->local_point_ids = NULL; /* Not needed for single-process */

  PLE_REALLOC(this_locator->distant_point_location,
              this_locator->n_interior,
              ple_lnum_t);
  PLE_REALLOC(this_locator->distant_point_coords,
              this_locator->n_interior * dim,
              ple_coord_t);

  n_interior = 0;
  n_exterior = 0;

  for (j = 0; j < n_points; j++) {

    if (point_list != NULL)
      coord_idx = point_list[j] - idb;
    else
      coord_idx = j;

    if (location[j] > -1) {
      this_locator->distant_point_location[n_interior] = location[j];
      for (l = 0; l < dim; l++) {
        this_locator->distant_point_coords[n_interior*dim + l]
          = point_coords[coord_idx*dim + l];
      }
      this_locator->interior_list[n_interior] = j + idb;
      n_interior += 1;
    }
    else {
      this_locator->exterior_list[n_exterior] = j + idb;
      n_exterior += 1;
    }

  }

}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *   interior      <-- if true, local_var is restricted to located points
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_distant(ple_locator_t     *this_locator,
                            void              *distant_var,
                            void              *local_var,
                            const ple_lnum_t  *local_list,
                            MPI_Datatype       datatype,
                            size_t             stride,
                            bool               reverse,
                            bool               interior)
{
  int dist_v_count, loc_v_count, size;
  int dist_rank;
  int dist_v_flag, loc_v_flag;
  ple_lnum_t n_points_loc, n_points_loc_max, n_points_dist;
  size_t dist_v_idx;
  void *dist_v_ptr;
  void *loc_v_buf;

  double comm_timing[4] = {0., 0., 0., 0.};

  MPI_Aint lb, extent;
  MPI_Status status;

  /* Check extent of datatype */

  MPI_Type_get_extent(datatype, &lb, &extent);
  MPI_Type_size(datatype, &size);

  if (extent != size)
    ple_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Initialization */

  n_points_loc_max = 0;

  for (int i = 0; i < this_locator->n_intersects; i++) {
    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];
    if (n_points_loc > n_points_loc_max)
      n_points_loc_max = n_points_loc;
  }

  PLE_MALLOC(loc_v_buf, n_points_loc_max*size*stride, char);

  const ple_lnum_t *il = this_locator->interior_list;

  /* Loop on MPI ranks */
  /*-------------------*/

  for (int li = 0; li < this_locator->n_intersects; li++) {

    int i = (this_locator->comm_order != NULL) ?
      this_locator->comm_order[li] : li;

    const ple_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (distant_var != NULL && n_points_dist > 0)
      dist_v_flag = 1;
    else
      dist_v_flag = 0;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Sendrecv(&dist_v_flag, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
                 &loc_v_flag, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
                 this_locator->comm, &status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    if (loc_v_flag == 1 && (local_var == NULL || n_points_loc == 0))
      ple_error(__FILE__, __LINE__, 0,
                _("Incoherent arguments to different instances in "
                  "_exchange_point_var().\n"
                  "Send and receive operations do not match "
                  "(dist_rank = %d\n)\n"), dist_rank);

    dist_v_idx = this_locator->distant_points_idx[i] * stride*size;
    dist_v_count = n_points_dist * stride * dist_v_flag;

    if (loc_v_flag > 0)
      loc_v_count = n_points_loc*stride;
    else
      loc_v_count = 0;

    /* Exchange information */

    if (distant_var != NULL)
      dist_v_ptr = (void *)(((char *)distant_var) + dist_v_idx);
    else
      dist_v_ptr = NULL;

    if (reverse == false) {

      _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   loc_v_buf, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          const size_t nbytes = stride*size;
          if (this_locator->n_exterior == 0 || interior) {
            for (int k = 0; k < n_points_loc; k++) {
              char *local_v_p = (char *)local_var + _local_point_ids[k]*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
          else {
            for (int k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                + il[_local_point_ids[k]]*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
        }
        else {
          if (this_locator->n_exterior == 0 || interior) {
            const size_t nbytes = stride*size;
            const ple_lnum_t idb = this_locator->point_id_base;
            for (int k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                + (local_list[_local_point_ids[k]] - idb)*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
          else {
            const size_t nbytes = stride*size;
            const ple_lnum_t idb = this_locator->point_id_base;
            for (int k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                +  (local_list[il[_local_point_ids[k]]] - idb)
                                  *nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
        }
      }

    }
    else { /* if (reverse == true) */

      if (loc_v_flag > 0) {
        if (local_list == NULL) {
          if (this_locator->n_exterior == 0 || interior) {
            const size_t nbytes = stride*size;
            for (int k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var + _local_point_ids[k]*nbytes;
              char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
          else {
            const size_t nbytes = stride*size;
            for (int k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var + il[_local_point_ids[k]]*nbytes;
              char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
        }
        else {
          const size_t nbytes = stride*size;
          const ple_lnum_t idb = this_locator->point_id_base;
          if (this_locator->n_exterior == 0 || interior) {
            for (int k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var
                  + (local_list[_local_point_ids[k]] - idb)*nbytes;
              char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
          else {
            for (int k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var
                  + (local_list[il[_local_point_ids[k]]] - idb)*nbytes;
              char *loc_v_buf_p = (char *)loc_v_buf + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
        }
      }

      _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

      MPI_Sendrecv(loc_v_buf, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                   this_locator->comm, &status);

      _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

    }

  } /* End of loop on MPI ranks */

  PLE_FREE(loc_v_buf);

  this_locator->exchange_wtime[1] += comm_timing[0];
  this_locator->exchange_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * This variant of the function uses asynchronous MPI calls.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list for local_var
 *   datatype      <-- variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *   interior      <-- if true, local_var is restricted to located points
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_distant_asyn(ple_locator_t     *this_locator,
                                 void              *distant_var,
                                 void              *local_var,
                                 const ple_lnum_t  *local_list,
                                 MPI_Datatype       datatype,
                                 size_t             stride,
                                 bool               reverse,
                                 bool               interior)
{
  int dist_v_count, loc_v_count, size;
  int dist_rank;
  ple_lnum_t n_points_loc, n_points_loc_tot, n_points_dist;
  size_t dist_v_idx;
  unsigned char *dist_v_ptr, *loc_v_ptr;

  MPI_Aint lb, extent;
  void *loc_v_buf = NULL;
  int *dist_v_flag = NULL, *loc_v_flag = NULL;
  MPI_Status *status = NULL;
  MPI_Request *request = NULL;

  double comm_timing[4] = {0., 0., 0., 0.};

  /* Check extent of datatype */

  MPI_Type_get_extent(datatype, &lb, &extent);
  MPI_Type_size(datatype, &size);

  if (extent != size)
    ple_error(__FILE__, __LINE__, 0,
              _("_exchange_point_var() is not implemented for use with\n"
                "MPI datatypes associated with structures using padding\n"
                "(for which size != extent)."));

  /* Initialization */

  n_points_loc_tot
    = this_locator->local_points_idx[this_locator->n_intersects];

  PLE_MALLOC(loc_v_flag, this_locator->n_intersects, int);
  PLE_MALLOC(dist_v_flag, this_locator->n_intersects, int);
  PLE_MALLOC(request, this_locator->n_intersects*2, MPI_Request);
  PLE_MALLOC(status, this_locator->n_intersects*2, MPI_Status);

  PLE_MALLOC(loc_v_buf, n_points_loc_tot*size*stride, char);

  const ple_lnum_t *il = this_locator->interior_list;

  /* First loop on distant ranks for argument checks */
  /*-------------------------------------------------*/

  _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

  for (int i = 0; i < this_locator->n_intersects; i++) {

    dist_rank = this_locator->intersect_rank[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (distant_var != NULL && n_points_dist > 0)
      dist_v_flag[i] = 1;
    else
      dist_v_flag[i] = 0;

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    MPI_Irecv(loc_v_flag + i, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
              this_locator->comm, &request[i*2]);
    MPI_Isend(dist_v_flag + i, 1, MPI_INT, dist_rank, PLE_MPI_TAG,
              this_locator->comm, &request[i*2+1]);
  }

  MPI_Waitall(this_locator->n_intersects*2, request, status);

  _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

  PLE_FREE(dist_v_flag);

  for (int i = 0; i < this_locator->n_intersects; i++) {

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =   this_locator->local_points_idx[i+1]
                   - this_locator->local_points_idx[i];

    if (loc_v_flag[i] == 1 && (local_var == NULL || n_points_loc == 0))
      ple_error(__FILE__, __LINE__, 0,
                _("Incoherent arguments to different instances in "
                  "_exchange_point_var().\n"
                  "Send and receive operations do not match "
                  "(dist_rank = %d\n)\n"), dist_rank);
  }

  /* Loop on distant ranks for exchange of data in standard mode */
  /*-------------------------------------------------------------*/

  if (reverse == false) {

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    loc_v_ptr = loc_v_buf;

    for (int i = 0; i < this_locator->n_intersects; i++) {

      dist_rank = this_locator->intersect_rank[i];

      n_points_loc =    this_locator->local_points_idx[i+1]
                      - this_locator->local_points_idx[i];

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;

      if (distant_var != NULL)
        dist_v_count = n_points_dist * stride;
      else
        dist_v_count = 0;

      if (loc_v_flag[i] > 0)
        loc_v_count = n_points_loc*stride;
      else
        loc_v_count = 0;

      /* Exchange information */

      if (distant_var != NULL)
        dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
      else
        dist_v_ptr = NULL;

      MPI_Irecv(loc_v_ptr, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2]);
      MPI_Isend(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2+1]);

      loc_v_ptr += loc_v_count*size;
    }

    MPI_Waitall(this_locator->n_intersects*2, request, status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);
  }

  /* Loop on distant ranks for preparation or retrieval of buffer data */
  /*-------------------------------------------------------------------*/

  loc_v_ptr = loc_v_buf;

  for (int i = 0; i < this_locator->n_intersects; i++) {

    const ple_lnum_t *_local_point_ids
      = this_locator->local_point_ids + this_locator->local_points_idx[i];

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
      - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    dist_v_idx = this_locator->distant_points_idx[i] * stride*size;

    if (distant_var != NULL)
      dist_v_count = n_points_dist * stride;
    else
      dist_v_count = 0;

    if (loc_v_flag[i] > 0)
      loc_v_count = n_points_loc*stride;
    else
      loc_v_count = 0;

    /* Exchange information */

    if (distant_var != NULL)
      dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
    else
      dist_v_ptr = NULL;

    dist_rank = this_locator->intersect_rank[i];

    n_points_loc =    this_locator->local_points_idx[i+1]
                    - this_locator->local_points_idx[i];

    n_points_dist =   this_locator->distant_points_idx[i+1]
                    - this_locator->distant_points_idx[i];

    if (reverse == false) {

      if (loc_v_flag[i] > 0) {
        if (local_list == NULL) {
          const size_t nbytes = stride*size;
          if (this_locator->n_exterior == 0 || interior) {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              char *local_v_p = (char *)local_var + _local_point_ids[k]*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
          else {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                + il[_local_point_ids[k]]*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
        }
        else {
          const size_t nbytes = stride*size;
          const ple_lnum_t idb = this_locator->point_id_base;
          if (this_locator->n_exterior == 0 || interior) {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                + (local_list[_local_point_ids[k]] - idb)*nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
          else {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              char *local_v_p =   (char *)local_var
                                +  (local_list[il[_local_point_ids[k]]] - idb)
                                  *nbytes;
              const char *loc_v_buf_p = (const char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                local_v_p[l] = loc_v_buf_p[l];
            }
          }
        }
      }

    }
    else { /* if (reverse == true) */

      if (loc_v_flag[i] > 0) {
        if (local_list == NULL) {
          const size_t nbytes = stride*size;
          if (this_locator->n_exterior == 0 || interior) {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var + _local_point_ids[k]*nbytes;
              char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
          else {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var + il[_local_point_ids[k]]*nbytes;
              char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
        }
        else {
          const size_t nbytes = stride*size;
          const ple_lnum_t idb = this_locator->point_id_base;
          if (this_locator->n_exterior == 0 || interior) {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var
                  + (local_list[_local_point_ids[k]] - idb)*nbytes;
              char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
          else {
            for (ple_lnum_t k = 0; k < n_points_loc; k++) {
              const char *local_v_p
                = (const char *)local_var
                  + (local_list[il[_local_point_ids[k]]] - idb)*nbytes;
              char *loc_v_buf_p = (char *)loc_v_ptr + k*nbytes;
              for (size_t l = 0; l < nbytes; l++)
                loc_v_buf_p[l] = local_v_p[l];
            }
          }
        }
      }

    }

    loc_v_ptr += loc_v_count*size;

  }

  /* Loop on distant ranks for exchange of data in reverse mode */
  /*------------------------------------------------------------*/

  if (reverse == true) {

    _locator_trace_start_comm(_ple_locator_log_start_p_comm, comm_timing);

    loc_v_ptr = loc_v_buf;

    for (int i = 0; i < this_locator->n_intersects; i++) {

      dist_rank = this_locator->intersect_rank[i];

      n_points_loc =    this_locator->local_points_idx[i+1]
                      - this_locator->local_points_idx[i];

      n_points_dist =   this_locator->distant_points_idx[i+1]
                      - this_locator->distant_points_idx[i];

      dist_v_idx = this_locator->distant_points_idx[i] * stride*size;

      if (distant_var != NULL)
        dist_v_count = n_points_dist * stride;
      else
        dist_v_count = 0;

      if (loc_v_flag[i] > 0)
        loc_v_count = n_points_loc*stride;
      else
        loc_v_count = 0;

      /* Exchange information */

      if (distant_var != NULL)
        dist_v_ptr = ((unsigned char *)distant_var) + dist_v_idx;
      else
        dist_v_ptr = NULL;

      MPI_Irecv(dist_v_ptr, dist_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2]);
      MPI_Isend(loc_v_ptr, loc_v_count, datatype, dist_rank, PLE_MPI_TAG,
                this_locator->comm, &request[i*2+1]);

      loc_v_ptr += loc_v_count*size;
    }

    MPI_Waitall(this_locator->n_intersects*2, request, status);

    _locator_trace_end_comm(_ple_locator_log_end_p_comm, comm_timing);

  }

  /* Free temporary arrays */

  PLE_FREE(loc_v_buf);

  PLE_FREE(loc_v_flag);
  PLE_FREE(request);
  PLE_FREE(status);

  this_locator->exchange_wtime[1] += comm_timing[0];
  this_locator->exchange_cpu_time[1] += comm_timing[1];
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Distribute variable defined on "distant points" to the original ("local")
 * points.
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_local(ple_locator_t     *this_locator,
                          void              *distant_var,
                          void              *local_var,
                          const ple_lnum_t  *local_list,
                          size_t             type_size,
                          size_t             stride,
                          bool               reverse)
{
  ple_lnum_t i;
  size_t j;
  ple_lnum_t n_points_loc;

  const size_t nbytes = stride*type_size;

  /* Initialization */

  if (this_locator->n_interior == 0)
    return;

  n_points_loc =   this_locator->local_points_idx[1]
                 - this_locator->local_points_idx[0];

  assert(n_points_loc == (  this_locator->distant_points_idx[1]
                          - this_locator->distant_points_idx[0]));

  /* Exchange information */

  if (reverse == false) {

    if (local_list == NULL)
      memcpy(local_var, distant_var, n_points_loc*nbytes);

    else {
      const ple_lnum_t idb = this_locator->point_id_base;
      for (i = 0; i < n_points_loc; i++) {
        char *local_var_p = (char *)local_var + (local_list[i] - idb)*nbytes;
        const char *distant_var_p = (const char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          local_var_p[j] = distant_var_p[j];
      }
    }

  }
  else { /* if (reverse == true) */

    if (local_list == NULL)
      memcpy(distant_var, local_var, n_points_loc*nbytes);

    else {
      const ple_lnum_t idb = this_locator->point_id_base;
      for (i = 0; i < n_points_loc; i++) {
        const char *local_var_p
          = (const char *)local_var + (local_list[i] - idb)*nbytes;
        char *distant_var_p = (char *)distant_var + i*nbytes;
        for (j = 0; j < nbytes; j++)
          distant_var_p[j] = local_var_p[j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Distribute variable defined on "distant points" to the original ("local")
 * points.
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *   local_var     <-> variable defined on local points (received)
 *   local_list    <-- optional indirection list for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if true, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

static void
_exchange_point_var_local_incomplete(ple_locator_t     *this_locator,
                                     void              *distant_var,
                                     void              *local_var,
                                     const ple_lnum_t  *local_list,
                                     size_t             type_size,
                                     size_t             stride,
                                     bool               reverse)
{
  const size_t nbytes = stride*type_size;

  const ple_lnum_t *il = this_locator->interior_list;

  /* Initialization */

  if (this_locator->n_interior == 0)
    return;

  ple_lnum_t n_points_loc =   this_locator->local_points_idx[1]
                            - this_locator->local_points_idx[0];

  assert(n_points_loc == (  this_locator->distant_points_idx[1]
                          - this_locator->distant_points_idx[0]));

  /* Exchange information */

  if (reverse == false) {

    if (local_list == NULL) {
      for (ple_lnum_t i = 0; i < n_points_loc; i++) {
        char *local_var_p = (char *)local_var + il[i]*nbytes;
        const char *distant_var_p = (const char *)distant_var + i*nbytes;
        for (size_t j = 0; j < nbytes; j++)
          local_var_p[j] = distant_var_p[j];
      }
    }
    else {
      const ple_lnum_t idb = this_locator->point_id_base;
      for (ple_lnum_t i = 0; i < n_points_loc; i++) {
        char *local_var_p = (char *)local_var + (local_list[il[i]] - idb)*nbytes;
        const char *distant_var_p = (const char *)distant_var + i*nbytes;
        for (size_t j = 0; j < nbytes; j++)
          local_var_p[j] = distant_var_p[j];
      }
    }

  }
  else { /* if (reverse == true) */

    if (local_list == NULL)
      for (ple_lnum_t i = 0; i < n_points_loc; i++) {
        const char *local_var_p = (const char *)local_var + il[i]*nbytes;
        char *distant_var_p = (char *)distant_var + i*nbytes;
        for (size_t j = 0; j < nbytes; j++)
          distant_var_p[j] = local_var_p[j];
      }

    else {
      const ple_lnum_t idb = this_locator->point_id_base;
      for (ple_lnum_t i = 0; i < n_points_loc; i++) {
        const char *local_var_p
          = (const char *)local_var + (local_list[il[i]] - idb)*nbytes;
        char *distant_var_p = (char *)distant_var + i*nbytes;
        for (size_t j = 0; j < nbytes; j++)
          distant_var_p[j] = local_var_p[j];
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute variable defined on distant points to processes owning
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
 * The local_var[] is defined at local points, restricted to located points
 * (those whose numbers are returned by ple_locator_get_interior_list())
 * if interior is set to 1.
 *
 * If the optional local_list indirection is used, it is assumed to use
 * the same base numbering as that defined by the options for the previous
 * call to ple_locator_set_mesh() or ple_locator_extend_search().
 *
 * \param[in]      this_locator pointer to locator structure
 * \param[in, out] distant_var  variable defined on distant points
 *                              (ready to send); size: n_dist_points*stride
 * \param[in, out] local_var    variable defined on local points
 *                              (received); size: n_interior*stride
 * \param[in]      local_list   optional indirection list for local_var
 * \param[in]      type_size    sizeof (float or double) variable type
 * \param[in]      stride       dimension (1 for scalar,
 *                              3 for interleaved vector)
 * \param[in]      reverse      if nonzero, exchange is reversed
 *                              (receive values associated with distant points
 *                              from the processes owning the original points)
 * \param[in]      interior     if true, local_var is restricted
 *                              to located points.
 */
/*----------------------------------------------------------------------------*/

static void
_exchange_point_var(ple_locator_t     *this_locator,
                    void              *distant_var,
                    void              *local_var,
                    const ple_lnum_t  *local_list,
                    size_t             type_size,
                    size_t             stride,
                    int                reverse,
                    bool               interior)
{
  double w_start, w_end, cpu_start, cpu_end;

  int mpi_flag = 0;
  bool _reverse = reverse;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

#if defined(PLE_HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    MPI_Datatype datatype = MPI_DATATYPE_NULL;

    if (type_size == sizeof(double))
      datatype = MPI_DOUBLE;
    else if (type_size == sizeof(float))
      datatype = MPI_FLOAT;
    else
      ple_error(__FILE__, __LINE__, 0,
                _("type_size passed to ple_locator_exchange_point_var() does\n"
                  "not correspond to double or float."));

    assert (datatype != MPI_DATATYPE_NULL);

    if (this_locator->exchange_algorithm == _EXCHANGE_SENDRECV)
      _exchange_point_var_distant(this_locator,
                                  distant_var,
                                  local_var,
                                  local_list,
                                  datatype,
                                  stride,
                                  _reverse,
                                  interior);

    else if (this_locator->exchange_algorithm == _EXCHANGE_ISEND_IRECV)
      _exchange_point_var_distant_asyn(this_locator,
                                       distant_var,
                                       local_var,
                                       local_list,
                                       datatype,
                                       stride,
                                       _reverse,
                                       interior);

  }

#endif /* defined(PLE_HAVE_MPI) */

  if (!mpi_flag) {
    if (this_locator->n_exterior == 0 || interior)
      _exchange_point_var_local(this_locator,
                                distant_var,
                                local_var,
                                local_list,
                                type_size,
                                stride,
                                _reverse);
    else
      _exchange_point_var_local_incomplete(this_locator,
                                           distant_var,
                                           local_var,
                                           local_list,
                                           type_size,
                                           stride,
                                           _reverse);
  }

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->exchange_wtime[0] += (w_end - w_start);
  this_locator->exchange_cpu_time[0] += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   time_type         <-- 0 for total times, 1 for communication times
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

static void
_get_times(const ple_locator_t  *this_locator,
           int                   time_type,
           double               *location_wtime,
           double               *location_cpu_time,
           double               *exchange_wtime,
           double               *exchange_cpu_time)
{
  const ple_locator_t  *_locator = this_locator;

  if (this_locator != NULL) {

    if (location_wtime != NULL) {
      *location_wtime = _locator->location_wtime[time_type];
    }
    if (location_cpu_time != NULL) {
      *location_cpu_time = _locator->location_cpu_time[time_type];
    }
    if (exchange_wtime != NULL)
      *exchange_wtime = _locator->exchange_wtime[time_type];
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = _locator->exchange_cpu_time[time_type];

  }
  else {

    if (location_wtime != NULL) {
      location_wtime[0] = 0.;
      location_wtime[1] = 0.;
    }
    if (location_cpu_time != NULL) {
      location_cpu_time[0] = 0.;
      location_cpu_time[1] = 0.;
    }
    if (exchange_wtime != NULL)
      *exchange_wtime = 0.;
    if (exchange_cpu_time != NULL)
      *exchange_cpu_time = 0.;
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Creation of a locator structure.
 *
 * Note that depending on the choice of ranks of the associated communicator,
 * distant ranks may in fact be truly distant or not. If n_ranks = 1 and
 * start_rank is equal to the current rank in the communicator, the locator
 * will work only locally.
 *
 * \param[in] comm       associated MPI communicator
 * \param[in] n_ranks    number of MPI ranks associated with distant location
 * \param[in] start_rank first MPI rank associated with distant location
 *
 * \return pointer to locator
 */
/*----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)
ple_locator_t *
ple_locator_create(MPI_Comm  comm,
                   int       n_ranks,
                   int       start_rank)
#else
ple_locator_t *
ple_locator_create(void)
#endif
{
  int  i;
  ple_locator_t  *this_locator;

  PLE_MALLOC(this_locator, 1, ple_locator_t);

  this_locator->dim = 0;
  this_locator->have_tags = 0;

#if defined(PLE_HAVE_MPI)
  this_locator->comm = comm;
  this_locator->n_ranks = n_ranks;
  this_locator->start_rank = start_rank;
#else
  this_locator->n_ranks = 1;
  this_locator->start_rank = 0;
#endif

  this_locator->locate_algorithm = _ple_locator_location_algorithm;
  this_locator->exchange_algorithm = _EXCHANGE_SENDRECV;
  this_locator->async_threshold = _ple_locator_async_threshold;

  this_locator->point_id_base = 0;

  this_locator->n_intersects = 0;
  this_locator->intersect_rank = NULL;
  this_locator->comm_order = NULL;

  this_locator->local_points_idx = NULL;
  this_locator->distant_points_idx = NULL;

  this_locator->local_point_ids = NULL;

  this_locator->distant_point_location = NULL;
  this_locator->distant_point_coords = NULL;

  this_locator->n_interior = 0;
  this_locator->interior_list = NULL;

  this_locator->n_exterior = 0;
  this_locator->exterior_list = NULL;

  for (i = 0; i < 2; i++) {
    this_locator->location_wtime[i] = 0.;
    this_locator->location_cpu_time[i] = 0.;
  }

  for (i = 0; i < 2; i++) {
    this_locator->exchange_wtime[i] = 0.;
    this_locator->exchange_cpu_time[i] = 0.;
  }

  return this_locator;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destruction of a locator structure.
 *
 * \param[in, out] this_locator locator to destroy
 *
 * \return NULL pointer
 */
/*----------------------------------------------------------------------------*/

ple_locator_t *
ple_locator_destroy(ple_locator_t  *this_locator)
{
  if (this_locator != NULL) {

    PLE_FREE(this_locator->local_points_idx);
    PLE_FREE(this_locator->distant_points_idx);

    if (this_locator->local_point_ids != NULL)
      PLE_FREE(this_locator->local_point_ids);

    PLE_FREE(this_locator->distant_point_location);
    PLE_FREE(this_locator->distant_point_coords);

    PLE_FREE(this_locator->intersect_rank);
    PLE_FREE(this_locator->comm_order);

    PLE_FREE(this_locator->interior_list);
    PLE_FREE(this_locator->exterior_list);

    PLE_FREE(this_locator);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare locator for use with a given mesh representation.
 *
 * \param[in, out] this_locator        pointer to locator structure
 * \param[in]      mesh                pointer to mesh representation structure
 * \param[in]      options             options array (size
 *                                     PLE_LOCATOR_N_OPTIONS), or NULL
 * \param[in]      tolerance_base      associated fixed tolerance
 * \param[in]      tolerance_fraction  associated fraction of element bounding
 *                                     boxes added to tolerance
 * \param[in]      dim                 spatial dimension of mesh and points to
 *                                     locate
 * \param[in]      n_points            number of points to locate
 * \param[in]      point_list          optional indirection array to point_coords
 * \param[in]      point_tag           optional point tag (size: n_points)
 * \param[in]      point_coords        coordinates of points to locate
 *                                     (dimension: dim * n_points)
 * \param[out]     distance            optional distance from point to matching
 *                                     element: < 0 if unlocated; 0 - 1 if inside
 *                                     and > 1 if outside a volume element, or
 *                                     absolute distance to a surface element
 *                                     (size: n_points)
 * \param[in]      mesh_extents_f      pointer to function computing mesh or mesh
 *                                     subset or element extents
 * \param[in]      mesh_locate_f       pointer to function wich updates the
 *                                     location[] and distance[] arrays
 *                                     associated with a set of points for
 *                                     points that are in an element of this
 *                                     mesh, or closer to one than to previously
 *                                     encountered elements.
 */
/*----------------------------------------------------------------------------*/

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
                     ple_mesh_elements_locate_t  *mesh_locate_f)
{
  double w_start, w_end, cpu_start, cpu_end;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

  /* Other initializations */

  this_locator->dim = dim;

  if (distance != NULL) {
    for (ple_lnum_t i = 0; i < n_points; i++)
      distance[i] = -1;
  }

  /* Release information if previously present */

  _clear_location_info(this_locator);

  ple_locator_extend_search(this_locator,
                            mesh,
                            options,
                            tolerance_base,
                            tolerance_fraction,
                            n_points,
                            point_list,
                            point_tag,
                            point_coords,
                            distance,
                            mesh_extents_f,
                            mesh_locate_f);

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->location_wtime[0] += (w_end - w_start);
  this_locator->location_cpu_time[0] += (cpu_end - cpu_start);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extend search for a locator for which set_mesh has already been
 *        called.
 *
 * \param[in, out] this_locator        pointer to locator structure
 * \param[in]      mesh                pointer to mesh representation structure
 * \param[in]      options             options array (size
 *                                     PLE_LOCATOR_N_OPTIONS), or NULL
 * \param[in]      tolerance_base      associated fixed tolerance
 * \param[in]      tolerance_fraction  associated fraction of element bounding
 *                                     boxes added to tolerance
 * \param[in]      n_points            number of points to locate
 * \param[in]      point_list          optional indirection array to point_coords
 * \param[in]      point_tag           optional point tag (size: n_points)
 * \param[in]      point_coords        coordinates of points to locate
 *                                     (dimension: dim * n_points)
 * \param[out]     distance            optional distance from point to matching
 *                                     element: < 0 if unlocated; 0 - 1 if inside
 *                                     and > 1 if outside a volume element, or
 *                                     absolute distance to a surface element
 *                                     (size: n_points)
 * \param[in]      mesh_extents_f      pointer to function computing mesh or mesh
 *                                     subset or element extents
 * \param[in]      mesh_locate_f       pointer to function wich updates the
 *                                     location[] and distance[] arrays
 *                                     associated with a set of points for
 *                                     points that are in an element of this
 *                                     mesh, or closer to one than to previously
 *                                     encountered elements.
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
                          const int                    point_tag[],
                          const ple_coord_t            point_coords[],
                          float                        distance[],
                          ple_mesh_extents_t          *mesh_extents_f,
                          ple_mesh_elements_locate_t  *mesh_locate_f)
{
  int i;
  double w_start, w_end, cpu_start, cpu_end;
  ple_lnum_t  *location;

  double comm_timing[4] = {0., 0., 0., 0.};
  int mpi_flag = 0;

  /* Initialize timing */

  w_start = ple_timer_wtime();
  cpu_start = ple_timer_cpu_time();

  if (options != NULL)
    this_locator->point_id_base = options[PLE_LOCATOR_NUMBERING];
  else
    this_locator->point_id_base = 0;

  const int idb = this_locator->point_id_base;

  this_locator->have_tags = 0;

  /* Prepare locator (MPI version) */
  /*-------------------------------*/

#if defined(PLE_HAVE_MPI)

  const int dim = this_locator->dim;

  MPI_Initialized(&mpi_flag);

  if (mpi_flag && this_locator->comm == MPI_COMM_NULL)
    mpi_flag = 0;

  if (mpi_flag) {

    /* Flag values
       0: mesh dimension
       1: space dimension
       2: minimum algorithm version
       3: maximum algorithm version
       4: preferred algorithm version
       5: have point tags */

    int globflag[6];
    int locflag[6] = {-1,
                      -1,
                      1, /* equivalent to _LOCATE_BB_SENDRECV */
                      -_LOCATE_BB_SENDRECV_ORDERED,
                      this_locator->locate_algorithm,
                      0};
    ple_lnum_t  *location_rank_id;

    /* Check that at least one of the local or distant nodal meshes
       is non-NULL, and at least one of the local or distant
       point sets is non null */

    if (mesh != NULL)
      locflag[0] = dim;

    if (n_points > 0)
      locflag[1] = dim;

    if (n_points > 0 && point_tag != NULL)
      locflag[5] = 1;

    _locator_trace_start_comm(_ple_locator_log_start_g_comm, comm_timing);

    MPI_Allreduce(locflag, globflag, 6, MPI_INT, MPI_MAX,
                  this_locator->comm);

    _locator_trace_end_comm(_ple_locator_log_end_g_comm, comm_timing);

    if (globflag[0] < 0 || globflag[1] < 0)
      return;
    else if (mesh != NULL && globflag[1] != dim)
      ple_error(__FILE__, __LINE__, 0,
                _("Locator trying to use distant space dimension %d\n"
                  "with local space dimension %d\n"),
                globflag[1], dim);
    else if (mesh == NULL && globflag[0] != dim)
      ple_error(__FILE__, __LINE__, 0,
                _("Locator trying to use local space dimension %d\n"
                  "with distant space dimension %d\n"),
                dim, globflag[0]);

    /* Check algorithm versions and supported features */

    globflag[3] = -globflag[3];

    /* Compatibility with older versions */
    for (i = 2; i < 5; i++) {
      if (globflag[i] == 1)
        globflag[i] = _LOCATE_BB_SENDRECV;
    }

    if (globflag[2] > globflag[3])
      ple_error(__FILE__, __LINE__, 0,
                _("Incompatible locator algorithm ranges:\n"
                  "  global minimum algorithm id %d\n"
                  "  global maximum algorithm id %d\n"
                  "PLE library versions or builds are incompatible."),
                globflag[2], globflag[3]);

    if (globflag[4] < globflag[2])
      globflag[4] = globflag[2];
    if (globflag[4] > globflag[3])
      globflag[4] = globflag[3];

    this_locator->locate_algorithm = globflag[4];

    if (globflag[5] > 0)
      this_locator->have_tags = 1;

    /* Free temporary memory */

    PLE_MALLOC(location, n_points, ple_lnum_t);
    PLE_MALLOC(location_rank_id, n_points, ple_lnum_t);

    _transfer_location_distant(this_locator,
                               n_points,
                               location,
                               location_rank_id);

    _locate_all_distant(this_locator,
                        mesh,
                        tolerance_base,
                        tolerance_fraction,
                        n_points,
                        point_list,
                        point_tag,
                        point_coords,
                        location,
                        location_rank_id,
                        distance,
                        mesh_extents_f,
                        mesh_locate_f);

    PLE_FREE(location_rank_id);
  }

#endif

  /* Prepare locator (local version) */
  /*---------------------------------*/

  if (!mpi_flag) {

    if (mesh == NULL || n_points == 0)
      return;

    if (point_tag != NULL)
      this_locator->have_tags = 1;

    PLE_MALLOC(location, n_points, ple_lnum_t);

    _transfer_location_local(this_locator,
                             n_points,
                             location);

    _locate_all_local(this_locator,
                      mesh,
                      tolerance_base,
                      tolerance_fraction,
                      n_points,
                      point_list,
                      point_tag,
                      point_coords,
                      location,
                      distance,
                      mesh_extents_f,
                      mesh_locate_f);

    PLE_FREE(location);

  }

  /* Update local_point_ids values */
  /*-------------------------------*/

  if (   this_locator->n_interior > 0
      && this_locator->local_point_ids != NULL) {

    ple_lnum_t  *reduced_index;

    PLE_MALLOC(reduced_index, n_points, ple_lnum_t);

    for (i = 0; i < n_points; i++)
      reduced_index[i] = -1;

    assert(  this_locator->local_points_idx[this_locator->n_intersects]
           == this_locator->n_interior);

    for (i = 0; i < this_locator->n_interior; i++)
      reduced_index[this_locator->interior_list[i] - idb] = i;

    /* Update this_locator->local_point_ids[] so that it refers
       to an index in a dense [0, this_locator->n_interior] subset
       of the local points */

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->local_point_ids[i]
        = reduced_index[this_locator->local_point_ids[i]];

    for (i = 0; i < this_locator->n_interior; i++)
      assert(this_locator->local_point_ids[i] > -1);

    PLE_FREE(reduced_index);

  }

  /* If an initial point list was given, update
     this_locator->interior_list and this_locator->exterior_list
     so that they refer to the same point set as that initial
     list (and not to an index within the selected point set) */

  if (point_list != NULL) {

    for (i = 0; i < this_locator->n_interior; i++)
      this_locator->interior_list[i]
        = point_list[this_locator->interior_list[i] - idb];

    for (i = 0; i < this_locator->n_exterior; i++)
      this_locator->exterior_list[i]
        = point_list[this_locator->exterior_list[i] - idb];

  }

  /* Finalize timing */

  w_end = ple_timer_wtime();
  cpu_end = ple_timer_cpu_time();

  this_locator->location_wtime[0] += (w_end - w_start);
  this_locator->location_cpu_time[0] += (cpu_end - cpu_start);

  this_locator->location_wtime[1] += comm_timing[0];
  this_locator->location_cpu_time[1] += comm_timing[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Shift location ids for located points after locator initialization.
 *
 * This is useful mainly to switch between 0-based to 1-based numberings.
 *
 * \param[in, out] this_locator    pointer to locator structure
 * \param[in]      location_shift  shift value
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_shift_locations(ple_locator_t  *this_locator,
                            ple_lnum_t      location_shift)
{
  int n_intersects = this_locator->n_intersects;
  if (n_intersects == 0)
    return;

  const ple_lnum_t n_points = this_locator->distant_points_idx[n_intersects];

  for (ple_lnum_t i = 0; i < n_points; i++) {
    if (this_locator->distant_point_location[i] > -1)
      this_locator->distant_point_location[i] += location_shift;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of distant points after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of distant points.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_dist_points(const ple_locator_t  *this_locator)
{
  ple_lnum_t retval = 0;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_points_idx[this_locator->n_intersects];
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an array of local element numbers containing (or nearest to)
 *  each distant point after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return local element numbers associated with distant points.
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_dist_locations(const ple_locator_t  *this_locator)
{
  const ple_lnum_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_ranks != 0)
      retval = this_locator->distant_point_location;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an array of coordinates of each distant point after
 * locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return coordinate array associated with distant points (interlaced).
 */
/*----------------------------------------------------------------------------*/

const ple_coord_t *
ple_locator_get_dist_coords(const ple_locator_t  *this_locator)
{
  const ple_coord_t * retval = NULL;

  if (this_locator != NULL) {
    if (this_locator->n_intersects != 0)
      retval = this_locator->distant_point_coords;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of points located after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of points located.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_interior(const ple_locator_t  *this_locator)
{
  ple_lnum_t retval = 0;

  if (this_locator != NULL)
    retval = this_locator->n_interior;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return list of points located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return list of points located.
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_interior_list(const ple_locator_t  *this_locator)
{
  return this_locator->interior_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of points not located after locator initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return number of points not located.
 */
/*----------------------------------------------------------------------------*/

ple_lnum_t
ple_locator_get_n_exterior(const ple_locator_t  *this_locator)
{
  return this_locator->n_exterior;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 *
 * \return list of points not located.
 */
/*----------------------------------------------------------------------------*/

const ple_lnum_t *
ple_locator_get_exterior_list(const ple_locator_t  *this_locator)
{
  return this_locator->exterior_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discard list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * \param[in] this_locator pointer to locator structure
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_discard_exterior(ple_locator_t  *this_locator)
{
  this_locator->n_exterior = 0;
  PLE_FREE(this_locator->exterior_list);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute variable defined on distant points to processes owning
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
 * \param[in]      this_locator pointer to locator structure
 * \param[in, out] distant_var  variable defined on distant points
 *                              (ready to send); size: n_dist_points*stride
 * \param[in, out] local_var    variable defined on located local points
 *                              (received); size: n_interior*stride
 * \param[in]      local_list   optional indirection list for local_var
 * \param[in]      type_size    sizeof (float or double) variable type
 * \param[in]      stride       dimension (1 for scalar,
 *                              3 for interleaved vector)
 * \param[in]      reverse      if nonzero, exchange is reversed
 *                              (receive values associated with distant points
 *                              from the processes owning the original points)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_exchange_point_var(ple_locator_t     *this_locator,
                               void              *distant_var,
                               void              *local_var,
                               const ple_lnum_t  *local_list,
                               size_t             type_size,
                               size_t             stride,
                               int                reverse)
{
  _exchange_point_var(this_locator,
                      distant_var,
                      local_var,
                      local_list,
                      type_size,
                      stride,
                      reverse,
                      true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute variable defined on distant points to processes owning
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
 * The local_var[] is defined at the local points (whether located or not)
 * provided when calling ple_locator_set_mesh() or ple_locator_extend_search().
 *
 * If the optional local_list indirection is used, it is assumed to use
 * the same base numbering as that defined by the options for the previous
 * call to ple_locator_set_mesh() or ple_locator_extend_search().
 *
 * \param[in]      this_locator pointer to locator structure
 * \param[in, out] distant_var  variable defined on distant points
 *                              (ready to send); size: n_dist_points*stride
 * \param[in, out] local_var    variable defined on local points
 *                              (received); size: n_points*stride
 * \param[in]      local_list   optional indirection list for local_var
 * \param[in]      type_size    sizeof (float or double) variable type
 * \param[in]      stride       dimension (1 for scalar,
 *                              3 for interleaved vector)
 * \param[in]      reverse      if nonzero, exchange is reversed
 *                              (receive values associated with distant points
 *                              from the processes owning the original points)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_exchange_point_var_all(ple_locator_t     *this_locator,
                                   void              *distant_var,
                                   void              *local_var,
                                   const ple_lnum_t  *local_list,
                                   size_t             type_size,
                                   size_t             stride,
                                   int                reverse)
{
  _exchange_point_var(this_locator,
                      distant_var,
                      local_var,
                      local_list,
                      type_size,
                      stride,
                      reverse,
                      false);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return timing information.
 *
 * In parallel mode, this includes communication time.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * \param[in]  this_locator      pointer to locator structure
 * \param[out] location_wtime    Location Wall-clock time (or NULL)
 * \param[out] location_cpu_time Location CPU time (or NULL)
 * \param[out] exchange_wtime    Variable exchange Wall-clock time
 *                               (size: 1 or NULL)
 * \param[out] exchange_cpu_time Variable exchange CPU time (size: 2 or NULL)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_get_times(const ple_locator_t  *this_locator,
                      double               *location_wtime,
                      double               *location_cpu_time,
                      double               *exchange_wtime,
                      double               *exchange_cpu_time)
{
  _get_times(this_locator,
             0,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return communication timing information.
 *
 * In serial mode, return times are always zero.
 *
 * When location on closest elements to force location of all points is
 * active, location times include a total value, followed by the value
 * associated with the location of closest elements stage.
 *
 * parameters:
 * \param[in]  this_locator      pointer to locator structure
 * \param[out] location_wtime    Location Wall-clock time (or NULL)
 * \param[out] location_cpu_time Location CPU time (or NULL)
 * \param[out] exchange_wtime    Variable exchange Wall-clock time (or NULL)
 * \param[out] exchange_cpu_time Variable exchange CPU time (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_get_comm_times(const ple_locator_t  *this_locator,
                           double               *location_wtime,
                           double               *location_cpu_time,
                           double               *exchange_wtime,
                           double               *exchange_cpu_time)
{
  _get_times(this_locator,
             1,
             location_wtime, location_cpu_time,
             exchange_wtime, exchange_cpu_time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump printout of a locator structure.
 *
 * \param this_locator pointer to structure that should be dumped
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_dump(const ple_locator_t  *this_locator)
{
  int  i;
  ple_lnum_t  j;
  const ple_lnum_t  *idx, *index, *loc;
  const ple_coord_t  *coords;

  const ple_locator_t  *_locator = this_locator;

  if (this_locator == NULL)
    return;

  /* Basic information */
  /*-------------------*/

  ple_printf("\n"
             "Locator:\n\n"
             "Spatial dimension:                     %d\n"
             "Exchange algorithm:                    %d\n"
             "Number of ranks of distant location:   %d\n"
             "First rank of distant location:        %d\n"
             "Number of intersecting distant ranks:  %d\n",
             _locator->dim,
             _locator->exchange_algorithm,
             _locator->n_ranks, _locator->start_rank,
             _locator->n_intersects);

#if defined(PLE_HAVE_MPI)
  if (_locator->comm != MPI_COMM_NULL)
    ple_printf("\n"
               "Associated MPI communicator:           %ld\n",
               (long)(_locator->comm));
#endif

  /* Arrays indexed by rank */
  /*------------------------*/

  if (_locator->intersect_rank != NULL) {
    for (i = 0; i < _locator->n_intersects; i++)
      ple_printf("\n"
                 "  Intersection %d with distant rank %d\n\n",
                 i, _locator->intersect_rank[i]);
  }
  if (_locator->intersect_rank != NULL) {
    for (i = 0; i < _locator->n_intersects; i++)
      ple_printf("\n"
                 "  Communication ordering %d: %d\n\n",
                 i, _locator->comm_order[i]);
  }

  if (_locator->n_interior > 0) {

    if (_locator->local_point_ids != NULL) {

      ple_printf("\n  Local point ids (for receiving):\n\n");
      idx = _locator->local_points_idx;
      index = _locator->local_point_ids;
      for (i = 0; i < _locator->n_intersects; i++) {
        if (idx[i+1] > idx[i]) {
          ple_printf("%6d (idx = %10d) %10d\n",
                     i, idx[i], index[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d\n", index[j]);
        }
        else {
          ple_printf("%6d (idx = %10d)\n", i, idx[i]);
        }
        ple_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
      }

    }

  }

  if (_locator->distant_points_idx != NULL) {

    idx = _locator->distant_points_idx;
    loc = _locator->distant_point_location;
    coords = _locator->distant_point_coords;

    if (idx[_locator->n_intersects] > 0)
      ple_printf("\n  Distant point location:\n\n");

    for (i = 0; i < _locator->n_intersects; i++) {

      if (idx[i+1] > idx[i]) {

        if (_locator->dim == 1) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e]\n",
                     i, _locator->intersect_rank[i], idx[i],
                     loc[idx[i]], coords[idx[i]]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d [%12.5e]\n",
                       loc[j], coords[j]);
        }
        else if (_locator->dim == 2) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e]\n",
                     i, idx[i], loc[idx[i]],
                     coords[2*idx[i]], coords[2*idx[i]+1]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          %10d [%12.5e, %12.5e]\n",
                       loc[j], coords[2*j], coords[2*j+1]);
        }
        else if (_locator->dim == 3) {
          ple_printf("%6d (idx = %10d) %10d [%12.5e, %12.5e, %12.5e]\n",
                     i, idx[i], loc[idx[i]],
                     coords[3*idx[i]], coords[3*idx[i]+1], coords[3*idx[i]+2]);
          for (j = idx[i] + 1; j < idx[i + 1]; j++)
            ple_printf("                          "
                       "%10d [%12.5e, %12.5e, %12.5e]\n",
                       loc[j], coords[3*j], coords[3*j+1], coords[3*j+2]);
        }

      } /* if (idx[i+1] > idx[i]) */

    }

    if (idx[_locator->n_intersects] > 0)
      ple_printf("   end (idx = %10d)\n", idx[_locator->n_intersects]);
  }

  /* Local arrays */
  /*--------------*/

  ple_printf("\n"
             "  Number of local points successfully located:  %d\n\n",
             _locator->n_interior);

  for (j = 0; j < _locator->n_interior; j++)
    ple_printf("    %10d\n", _locator->interior_list[j]);

  if  (_locator->n_interior > 0)
    ple_printf("\n");

  ple_printf("  Number of local points not located:  %d\n",
             _locator->n_exterior);

  for (j = 0; j < _locator->n_exterior; j++)
    ple_printf("    %10d\n", _locator->exterior_list[j]);

  if  (_locator->n_exterior > 0)
    ple_printf("\n");

  /* Timing information */
  /*--------------------*/

  ple_printf("  Location Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->location_wtime[0], _locator->location_wtime[1]);

  ple_printf("  Location CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->location_cpu_time[0], _locator->location_cpu_time[1]);

  ple_printf("  Exchange Wall-clock time: %12.5f (comm: %12.5f)\n",
             _locator->exchange_wtime[0], _locator->exchange_wtime[1]);

  ple_printf("  Exchange CPU time:        %12.5f (comm: %12.5f)\n",
             _locator->exchange_cpu_time[0], _locator->exchange_cpu_time[1]);

}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the maximum number of exchanging ranks for which we use
 * asynchronous MPI sends and receives instead of MPI_SendRecv.
 *
 * \return the maximum number of ranks allowing asynchronous exchanges
 */
/*----------------------------------------------------------------------------*/

int
ple_locator_get_async_threshold(void)
{
  return _ple_locator_async_threshold;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the maximum number of exchanging ranks for which we use
 * asynchronous MPI sends and receives instead of MPI_SendRecv.
 *
 * \param threshold maximum number of ranks allowing asynchronous exchanges
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_async_threshold(int threshold)
{
  _ple_locator_async_threshold = threshold;
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default locator options.
 *
 * \param[in]  key    option name
 * \param[in]  value  associated value
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_default_option(const char  *key,
                               const char  *value)
{
  _parse_locator_option(key,
                        value,
                        &_ple_locator_location_algorithm,
                        &_ple_locator_async_threshold);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set options for a given locator.
 *
 * \param[in, out]  this_locator  pointer to locator structure
 * \param[in]       key           option name
 * \param[in]       value         associated value
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_options(ple_locator_t  *this_locator,
                        const char     *key,
                        const char     *value)
{
  if (key == NULL || value == NULL)
    return;

  _parse_locator_option(key,
                        value,
                        &(this_locator->locate_algorithm),
                        &(this_locator->async_threshold));
}

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Register communication logging functions for locator instrumentation.
 *
 * By default, locators are not instrumented.
 *
 * Functions using MPE may be defined and used, but other similar systems
 * may be used.
 *
 * \param[in] log_function pointer to logging function
 * \param[in] start_p_comm point to point communication start event number
 * \param[in] end_p_comm   point to point communication end event number
 * \param[in] start_g_comm global communication start event number
 * \param[in] end_g_comm   global communication end event number
 */
/*----------------------------------------------------------------------------*/

void
ple_locator_set_comm_log(ple_locator_log_t  *log_function,
                         int                 start_p_comm,
                         int                 end_p_comm,
                         int                 start_g_comm,
                         int                 end_g_comm)
{
  _ple_locator_log_func = log_function;

  _ple_locator_log_start_p_comm = start_p_comm;
  _ple_locator_log_end_p_comm = end_p_comm;
  _ple_locator_log_start_g_comm = start_g_comm;
  _ple_locator_log_end_g_comm = end_g_comm;
}

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
