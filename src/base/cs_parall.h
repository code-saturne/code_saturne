#ifndef __CS_PARALL_H__
#define __CS_PARALL_H__

/*============================================================================
 * Functions dealing with parallelism
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include "cs_execution_context.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * General types and macros used throughout code_saturne
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Variable value type.
 *----------------------------------------------------------------------------*/

/*! Algorithm classes for indirect sums from a graph's edges to its nodes,
  such as face-to cell sums. */

typedef enum {

  CS_E2N_SUM_SCATTER,           /*!< Iterate on edges and scatter values to
                                  nodes, assuming no data races
                                  (serial / single thread mode). */
  CS_E2N_SUM_SCATTER_ATOMIC,    /*!< Iterate on edges and scatter values to
                                  nodes, using atomic sums. */
  CS_E2N_SUM_GATHER             /*!< Iterate on nodes and gather (indexed)
                                  edge values */

} cs_e2n_sum_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Preferred indexed sum option, adapted to shared-memory parallelism */

extern cs_e2n_sum_t cs_glob_e2n_sum_type;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a counter on all default communicator processes.
 *
 * \param[in, out]  cpt  local counter in, global counter out (size: n)
 * \param[in]       n    number of values
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_counter(cs_gnum_t   cpt[],
                  const int   n)
{
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, cpt, n, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
}

#else

void
cs_parall_counter(cs_gnum_t   cpt[],
                  const int   n);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a counter on all default communicator processes.
 *
 * \param[in, out]  cpt  local counter in, global counter out (size: n)
 * \param[in]       n    number of values
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_counter_max(cs_lnum_t   cpt[],
                      const int   n)
{
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, cpt, n, CS_MPI_LNUM, MPI_MAX,
                  cs_glob_mpi_comm);
  }
}

#else

void
cs_parall_counter_max(cs_lnum_t   cpt[],
                      const int   n);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a given datatype on all default communicator processes.
 *
 * \param[in]       n        number of values
 * \param[in]       datatype matching code_saturne datatype
 * \param[in, out]  val      local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_sum(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, val, n, cs_datatype_to_mpi[datatype], MPI_SUM,
                  cs_glob_mpi_comm);
  }
}

#else

void
cs_parall_sum(int             n,
              cs_datatype_t   datatype,
              void           *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype on all default
 *        communicator processes.
 *
 * \param[in]       n        number of values
 * \param[in]       datatype matching code_saturne datatype
 * \param[in, out]  val      local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_max(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, val, n, cs_datatype_to_mpi[datatype], MPI_MAX,
                  cs_glob_mpi_comm);
  }
}

#else

void
cs_parall_max(int             n,
              cs_datatype_t   datatype,
              void           *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype on all default
 *        communicator processes.
 *
 * \param[in]       n        number of values
 * \param[in]       datatype matching code_saturne datatype
 * \param[in, out]  val      local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_min(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, val, n, cs_datatype_to_mpi[datatype], MPI_MIN,
                  cs_glob_mpi_comm);
  }
}

#else

void
cs_parall_min(int             n,
              cs_datatype_t   datatype,
              void           *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Broadcast values of a given datatype to all
 *        default communicator processes.
 *
 * \param[in]       root_rank rank from which to broadcast
 * \param[in]       n         number of values
 * \param[in]       datatype  matching code_saturne datatype
 * \param[in, out]  val       values to broadcast; input on root_rank,
 *                            output on others (size: n)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

inline static void
cs_parall_bcast(int             root_rank,
                int             n,
                cs_datatype_t   datatype,
                void           *val)
{
  if (cs_glob_n_ranks > 1)
    MPI_Bcast(val, n, cs_datatype_to_mpi[datatype], root_rank,
              cs_glob_mpi_comm);
}

#else

#define cs_parall_bcast(_root_rank, _n, _datatype, _val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a global array from each local array in each domain.
 *
 * Local arrays are appended in order of owning MPI rank.
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   n_elts    size of the local array
 * \param[in]   n_g_elts  size of the global array
 * \param[in]   array     local array (size: n_elts)
 * \param[out]  g_array   global array  (size: n_g_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_allgather_r(int        n_elts,
                      int        n_g_elts,
                      cs_real_t  array[],
                      cs_real_t  g_array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build an ordered global array from each local array in each domain.
 *
 * Local array elements are ordered based on a given key value (usually
 * some form of coordinate, so the result should be independent of
 * partitioning as long as the definition of the o_key array's defintion
 * is itself independent of the partitioning.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   n_elts    number of local elements
 * \param[in]   n_g_elts  number of global elements
 * \param[in]   stride    number of values per element
 * \param[in]   o_key     ordering key (coordinate) value per element
 * \param[in]   array     local array (size: n_elts*stride)
 * \param[out]  g_array   global array  (size: n_g_elts*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_allgather_ordered_r(int        n_elts,
                              int        n_g_elts,
                              int        stride,
                              cs_real_t  o_key[],
                              cs_real_t  array[],
                              cs_real_t  g_array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a global array on the given root rank from all local arrays.
 *
 * Local arrays are appended in order of owning MPI rank.
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     size of the local array
 * \param[in]   n_g_elts   size of the global array
 * \param[in]   array      local array (size: n_elts)
 * \param[out]  g_array    global array  (size: n_g_elts) only usable by the
 *                         root rank
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_gather_r(int               root_rank,
                   int               n_elts,
                   int               n_g_elts,
                   const cs_real_t   array[],
                   cs_real_t         g_array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build an ordered global array on the given root rank from
 *        all local arrays.
 *
 * Local array elements are ordered based on a given key value (usually
 * some form of coordinate, so the result should be independent of
 * partitioning as long as the definition of the o_key array's defintion
 * is itself independent of the partitioning.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     number of local elements
 * \param[in]   n_g_elts   number of global elements
 * \param[in]   stride     number of values per element
 * \param[in]   o_key      ordering key (coordinate) value per element
 * \param[in]   array      local array (size: n_elts*stride)
 * \param[out]  g_array    global array  (size: n_g_elts*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_gather_ordered_r(int        root_rank,
                           int        n_elts,
                           int        n_g_elts,
                           int        stride,
                           cs_real_t  o_key[],
                           cs_real_t  array[],
                           cs_real_t  g_array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute a global array from a given root rank over all ranks.
 *        Each rank receive the part related to its local elements.
 *
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to specific usage, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     size of the local array
 * \param[in]   n_g_elts   size of the global array
 * \param[in]   g_array    global array  (size: n_g_elts) only usable by the
 *                         root rank
 * \param[out]  array      local array (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_scatter_r(int               root_rank,
                    int               n_elts,
                    int               n_g_elts,
                    const cs_real_t   g_array[],
                    cs_real_t         array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a global array on the given root rank from all local arrays.
 *        Function dealing with single-precision arrays.
 *
 * Local arrays are appended in order of owning MPI rank.
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     size of the local array
 * \param[in]   n_g_elts   size of the global array
 * \param[in]   array      local array (size: n_elts)
 * \param[out]  g_array    global array  (size: n_g_elts) only usable by the
 *                         root rank
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_gather_f(int             root_rank,
                   int             n_elts,
                   int             n_g_elts,
                   const float     array[],
                   float           g_array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute a global array from a given root rank over all ranks.
 *        Each rank receive the part related to its local elements.
 *        Function dealing with single-precision arrays.
 *
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to specific usage, as it may limit scalability (especially as regards
 * memory usage).
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     size of the local array
 * \param[in]   n_g_elts   size of the global array
 * \param[in]   g_array    global array  (size: n_g_elts) only usable by the
 *                         root rank
 * \param[out]  array      local array (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_scatter_f(int           root_rank,
                    int           n_elts,
                    int           n_g_elts,
                    const float   g_array[],
                    float         array[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum value of a real and the value of related array on all
 *        default communicator processes.
 *
 * \param[in]       n             size of the related array
 * \param[in, out]  max           local max in, global max out
 * \param[in, out]  max_loc_vals  array values at location of local max in,
 *                                and at location of global max out
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_max_loc_vals(int         n,
                       cs_real_t  *max,
                       cs_real_t   max_loc_vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum value of a real and the value of related array on all
 *        default communicator processes.
 *
 * \param[in]       n             size of the related array
 * \param[in, out]  minx           local max in, global max out
 * \param[in, out]  min_loc_vals  array values at location of local max in,
 *                                and at location of global max out
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_min_loc_vals(int         n,
                       cs_real_t  *min,
                       cs_real_t   min_loc_vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Given an (id, rank, value) tuple, return the local id and rank
 *        corresponding to the global minimum value.
 *
 * \param[in, out]   elt_id   element id for which the value is the smallest
 *                            (local in, global out)
 * \param[in, out]   rank_id  rank id for which the value is the smallest
 *                            (local in, global out)
 * \param[in]        val      associated local minimum value
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_min_id_rank_r(cs_lnum_t  *elt_id,
                        int        *rank_id,
                        cs_real_t   val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return minimum recommended scatter or gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * \return  minimum recommended part to block or gather buffer size (in bytes)
 */
/*----------------------------------------------------------------------------*/

size_t
cs_parall_get_min_coll_buf_size(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define minimum recommended scatter or gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * \param[in]  buffer_size  minimum recommended part to block
 *             or gather buffer size (in bytes)
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_set_min_coll_buf_size(size_t buffer_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute recommended number of threads for a section.
 *
 * This is based on the available number of threads, and limited
 * by a minimum number of elements per thread.
 *
 * \param[in]  n_elements           size of array
 * \param[in]  min_thread_elements  minimum number of elements per thread.
 */
/*----------------------------------------------------------------------------*/

inline static int
cs_parall_n_threads(cs_lnum_t  n_elements,
                    cs_lnum_t  min_thread_elements)
{
#if defined(HAVE_OPENMP)
  int n_t = omp_get_max_threads();
  int n_t_l = n_elements / min_thread_elements;
  if (n_t_l < n_t)
    n_t = n_t_l;
  if (n_t < 1)
    n_t = 1;
  return n_t;
#else
  CS_UNUSED(n_elements);         /* avoid compiler warning */
  CS_UNUSED(min_thread_elements);
  return 1;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *        When called inside an OpenMP parallel section, this will return the
 *        start an past-the-end indexes for the array range assigned to that
 *        thread. In other cases, the start index is 1, and the past-the-end
 *        index is n;
 *
 * \param[in]       n          size of array
 * \param[in]       type_size  element type size (or multiple)
 * \param[in, out]  s_id       start index for the current thread
 * \param[in, out]  e_id       past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_parall_thread_range(cs_lnum_t    n,
                       size_t       type_size,
                       cs_lnum_t   *s_id,
                       cs_lnum_t   *e_id)
{
#if defined(HAVE_OPENMP)
  const int t_id = omp_get_thread_num();
  const int n_t = omp_get_num_threads();
  const cs_lnum_t t_n = (n + n_t - 1) / n_t;
  const cs_lnum_t cl_m = CS_CL_SIZE / type_size;  /* Cache line multiple */

  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, cl_m);
  *e_id = cs_align(*e_id, cl_m);
  if (*e_id > n) *e_id = n;
#else
  CS_UNUSED(type_size);         /* avoid compiler warning */
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread for upper triangular
 *        matrix elements.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that
 * thread. In other cases, the start index is 1, and the past-the-end
 * index is n;
 *
 * Compared to \ref cs_parall_thread_range, this variant assumes work on
 * the upper triangular part of a matrix, where the lower part is
 * ignored.
 *
 * \param[in]       n          size of array
 * \param[in]       type_size  element type size (or multiple)
 * \param[in, out]  s_id       start index for the current thread
 * \param[in, out]  e_id       past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_parall_thread_range_upper(cs_lnum_t    n,
                             size_t       type_size,
                             cs_lnum_t   *s_id,
                             cs_lnum_t   *e_id)
{
#if defined(HAVE_OPENMP)
  const int t_id = omp_get_thread_num();
  const double n_t = omp_get_num_threads();
  const cs_lnum_t cl_m = CS_CL_SIZE / type_size;  /* Cache line multiple */

  double r0 = (double)t_id / (double)n_t;
  double r1 = (double)(t_id+1) / (double)n_t;

  r0 = r0*r0;
  r1 = r1*r1;

  const cs_lnum_t t_0 = r0*n;
  const cs_lnum_t t_1 = r1*n;

  *s_id = t_0 * n;
  *e_id = t_1 * n;
  *s_id = cs_align(*s_id, cl_m);
  *e_id = cs_align(*e_id, cl_m);
  if (*e_id > n) *e_id = n;
#else
  CS_UNUSED(type_size);         /* avoid compiler warning */
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute number of blocks needed for a given array and block sizes.
 *
 * \param[in]  n           size of array
 * \param[in]  block_size  block size for sub-loops
 *
 * \return  number of required blocks
 */
/*----------------------------------------------------------------------------*/

static inline size_t
cs_parall_block_count(size_t  n,
                      size_t  block_size)
{
  return (n % block_size) ?  n/block_size + 1 : n/block_size;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#if defined(__cplusplus)

/*=============================================================================
 * Public C++ functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a counter on all default communicator processes.
 *
 * \param[in]       ec   parallel execution context
 * \param[in, out]  cpt  local counter in, global counter out (size: n)
 * \param[in]       n    number of values
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_counter(cs_execution_context  *ec,
                  cs_gnum_t              cpt[],
                  const int              n)
{
  if (ec->use_mpi()) {
    MPI_Allreduce(MPI_IN_PLACE, cpt, n, CS_MPI_GNUM, MPI_SUM,
                  ec->comm());
  }
}

#else

void
cs_parall_counter(cs_execution_context  *ec,
                  cs_gnum_t              cpt[],
                  const int              n);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a counter on all default communicator processes.
 *
 * \param[in]       ec   parallel execution context
 * \param[in, out]  cpt  local counter in, global counter out (size: n)
 * \param[in]       n    number of values
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_counter_max(cs_execution_context  *ec,
                      cs_lnum_t              cpt[],
                      const int              n)
{
  if (ec->use_mpi()) {
    MPI_Allreduce(MPI_IN_PLACE, cpt, n, CS_MPI_LNUM, MPI_MAX,
                  ec->comm());
  }
}

#else

void
cs_parall_counter_max(cs_execution_context  *ec,
                      cs_lnum_t              cpt[],
                      const int              n);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a given datatype on all of a communicator processes.
 *
 * \param[in]       ec        parallel execution context
 * \param[in]       n         number of values
 * \param[in]       datatype  matching code_saturne datatype
 * \param[in, out]  val       local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_sum(cs_execution_context  *ec,
              int                    n,
              cs_datatype_t          datatype,
              void                  *val)
{
  if (ec->use_mpi()) {
    MPI_Allreduce(MPI_IN_PLACE, val, n,
                  cs_datatype_to_mpi[datatype], MPI_SUM,
                  ec->comm());
  }
}

#else

void
cs_parall_sum(cs_execution_context  *ec,
              int                    n,
              cs_datatype_t          datatype,
              void                  *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype on all default
 *        communicator processes.
 *
 * \param[in]       ec       parallel execution context
 * \param[in]       n        number of values
 * \param[in]       datatype matching code_saturne datatype
 * \param[in, out]  val      local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_max(cs_execution_context  *ec,
              int                    n,
              cs_datatype_t          datatype,
              void                  *val)
{
  if (ec->use_mpi()) {
    MPI_Allreduce(MPI_IN_PLACE, val, n,
                  cs_datatype_to_mpi[datatype], MPI_MAX,
                  ec->comm());
  }
}

#else

void
cs_parall_max(cs_execution_context *ec,
              int                   n,
              cs_datatype_t         datatype,
              void                 *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype on all default
 *        communicator processes.
 *
 * \param[in]       ec       parallel execution context
 * \param[in]       n        number of values
 * \param[in]       datatype matching code_saturne datatype
 * \param[in, out]  val      local value input, global value output (array)
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI_IN_PLACE)

inline static void
cs_parall_min(cs_execution_context  *ec,
              int                    n,
              cs_datatype_t          datatype,
              void                  *val)
{
  if (ec->use_mpi()) {
    MPI_Allreduce(MPI_IN_PLACE, val, n,
                  cs_datatype_to_mpi[datatype], MPI_MIN,
                  ec->comm());
  }
}

#else

void
cs_parall_min(cs_execution_context  *ec,
              int                    n,
              cs_datatype_t          datatype,
              void                  *val);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a given thread.
 *        This will return the
 *        start an past-the-end indexes for the array range assigned to that
 *        thread. In other cases, the start index is 1, and the past-the-end
 *        index is n;
 *
 * \param[in]       n          size of array
 * \param[in]       type_size  element type size (or multiple)
 * \param[in]       t_id       thread id
 * \param[in]       n_t        number of threads
 * \param[in, out]  s_id       start index for the current thread
 * \param[in, out]  e_id       past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_parall_thread_range(cs_lnum_t    n,
                       size_t       type_size,
                       int          t_id,
                       int          n_t,
                       cs_lnum_t   *s_id,
                       cs_lnum_t   *e_id)
{
  const cs_lnum_t t_n = (n + n_t - 1) / n_t;
  const cs_lnum_t cl_m = CS_CL_SIZE / type_size;  /* Cache line multiple */

  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, cl_m);
  *e_id = cs_align(*e_id, cl_m);
  if (*e_id > n) *e_id = n;
}

/*=============================================================================
 * Public C++ templates
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a given datatype on all default communicator processes.
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_sum_scalars
(
  T&       first, /*!< [in, out]  First scalar to update */
  Vals&... values /*!< [in, out]  Additional scalars to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_sum(1, datatype, &first);
  else {
    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i+1] = *(_values[i]);

    cs_parall_sum(n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i+1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum values of a given datatype over a given communicator
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_sum_scalars
(
  cs_execution_context  *ec,     /*!<[in] Parallel execution context */
  T&                     first,  /*!<[in, out]  First scalar to update */
  Vals&...               values  /*!<[in, out]  Additional scalars to update */
)
{
#if defined(HAVE_MPI)

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_sum(ec, 1, datatype, &first);
  else {
    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i+1] = *(_values[i]);

    cs_parall_sum(ec, n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i+1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum strided-values of a given datatype over a communicator
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

template <int Stride, typename T, typename... Vals>
static void
cs_parall_sum_strided
(
  cs_execution_context  *ec,      /*!< [in] Parallel execution context */
  T                      first[], /*!< [in, out]  First scalar to update */
  Vals&&...              values   /*!< [in, out]  Additional values to update */
)
{
#if defined(HAVE_MPI)

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_sum(ec, Stride, datatype, first);
  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_sum(ec, work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++) {
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
    }
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sum strided-values of a given datatype on all default communicator
 *        processes.
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

template <int Stride, typename T, typename... Vals>
static void
cs_parall_sum_strided
(
  T         first[], /*!< [in, out]  First scalar to update */
  Vals&&... values   /*!< [in, out]  Additional values to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_sum(Stride, datatype, first);
  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_sum(work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++) {
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
    }
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype on all default
 *        communicator processes.
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_max_scalars
(
  T&       first, /*!< [in, out]  First scalar to update */
  Vals&... values /*!< [in, out]  Additional scalars to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_max(1, datatype, &first);
  else {

    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i+1] = *(_values[i]);

    cs_parall_max(n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i+1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype over a communicator.
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_max_scalars
(
  cs_execution_context  *ec,    /*!< [in] Parallel execution context */
  T&                     first, /*!< [in, out]  First scalar to update */
  Vals&...               values /*!< [in, out]  Additional scalars to update */
)
{
#if defined(HAVE_MPI)

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_max(ec, 1, datatype, &first);
  else {

    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i+1] = *(_values[i]);

    cs_parall_max(ec, n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i+1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype on all default
 *        communicator processes.
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

template <int Stride, typename T, typename... Vals>
static void
cs_parall_max_strided
(
  T         first[], /*!< [in, out]  First scalar to update */
  Vals&&... values   /*!< [in, out]  Additional values to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_max(Stride, datatype, first);
  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_max(work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum values of a given datatype over a communicator.
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

template <int Stride, typename T, typename... Vals>
static void
cs_parall_max_strided
(
  cs_execution_context  *ec,      /*!< [in] Parallel execution context */
  T                      first[], /*!< [in, out]  First scalar to update */
  Vals&&...              values   /*!< [in, out]  Additional values to update */
)
{
#if defined(HAVE_MPI)

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_max(ec, Stride, datatype, first);
  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_max(ec, work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype on all default
 *        communicator processes.
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_min_scalars
(
  T&       first, /*!< [in, out]  First scalar to update */
  Vals&... values /*!< [in, out]  Additional scalars to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  if (n_vals == 0)
    cs_parall_min(1, datatype, &first);

  else {
    /* Temporary work array and parallel sum */

    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i + 1] = *(_values[i]);

    cs_parall_min(n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i + 1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype over a communicator.
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename... Vals>
static void
cs_parall_min_scalars
(
  cs_execution_context  *ec,     /*!< [in] Parallel execution context */
  T&                     first,  /*!< [in, out]  First scalar to update */
  Vals&...               values  /*!< [in, out]  Additional scalar to update */
)
{
#if defined(HAVE_MPI)

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  if (n_vals == 0)
    cs_parall_min(ec, 1, datatype, &first);

  else {
    /* Temporary work array and parallel sum */

    /* Unpack values */
    T *_values[] = {&values ...};

    T w[n_vals + 1];
    w[0] = first;
    for (size_t i = 0; i < n_vals; i++)
      w[i + 1] = *(_values[i]);

    cs_parall_min(ec, n_vals + 1, datatype, w);

    first = w[0];
    for (size_t i = 0; i < n_vals; i++)
      *(_values[i]) = w[i + 1];
  }

#endif // defined(HAVE_MPI)
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype on all default
 *        communicator processes.
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

template <int Stride, typename T, typename... Vals>
static void
cs_parall_min_strided
(
  T         first[],  /*!< [in, out]  First value to update */
  Vals&&... values    /*!< [in, out]  Additional values to update */
)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks == 1)
    return;

  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_min(Stride, datatype, first);
  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_min(work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum values of a given datatype over a communicator.
 *
 * \tparam T      : datatype
 * \tparam stride : stride/dimension of values
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

template <int Stride, typename T, typename... Vals>
static void
cs_parall_min_strided
(
  cs_execution_context  *ec,      /*!< [in] Parallel execution context */
  T                      first[], /*!< [in, out]  First value to update */
  Vals&&...              values   /*!< [in, out]  Additional values to update */
)
{
  /* Count number of values */
  constexpr size_t n_vals = sizeof...(Vals);

  /* Set datatype for global communication */
  cs_datatype_t datatype = cs_datatype_from_type<T>();

  /* Temporary work array and parallel sum */
  if (n_vals == 0)
    cs_parall_min(ec, Stride, datatype, first);

  else {
    /* Unpack values */
    T *_values[] = {values ...};

    constexpr size_t work_size = (n_vals + 1) * Stride;

    T w[work_size];
    for (int i = 0; i < Stride; i++)
      w[i] = first[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        w[(i+1)*Stride + j] = _values[i][j];

    cs_parall_min(ec, work_size, datatype, w);

    for (int i = 0; i < Stride; i++)
      first[i] = w[i];

    for (size_t i = 0; i < n_vals; i++)
      for (int j = 0; j < Stride; j++)
        _values[i][j] = w[(i+1)*Stride + j];
  }
}

#endif // defined(HAVE_MPI)

#endif //__cplusplus

/*----------------------------------------------------------------------------*/

#endif /* __CS_PARALL_H__ */
