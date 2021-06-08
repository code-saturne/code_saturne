#ifndef __CS_PARALL_H__
#define __CS_PARALL_H__

/*============================================================================
 * Functions dealing with parallelism
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Sum values of a counter on all default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter in, global counter out (size: n)
 *   n   <-- number of values
 *----------------------------------------------------------------------------*/

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

#elif defined(HAVE_MPI)

void
cs_parall_counter(cs_gnum_t   cpt[],
                  const int   n);

#else

#define cs_parall_counter(_cpt, _n)

#endif

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter in, global counter out (size: n)
 *   n   <-> number of values
 *----------------------------------------------------------------------------*/

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

#elif defined(HAVE_MPI)

void
cs_parall_counter_max(cs_lnum_t   cpt[],
                      const int   n);

#else

#define cs_parall_counter_max(_cpt, _n)

#endif

/*----------------------------------------------------------------------------
 * Sum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local sum in, global sum out (array)
 *----------------------------------------------------------------------------*/

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

#elif defined(HAVE_MPI)

void
cs_parall_sum(int             n,
              cs_datatype_t   datatype,
              void           *val);

#else

#define cs_parall_sum(_n, _datatype, _val) { };

#endif

/*----------------------------------------------------------------------------
 * Maximum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

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

#elif defined(HAVE_MPI)

void
cs_parall_max(int             n,
              cs_datatype_t   datatype,
              void           *val);

#else

#define cs_parall_max(_n, _datatype, _val);

#endif

/*----------------------------------------------------------------------------
 * Minimum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching Code_Saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

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

#elif defined(HAVE_MPI)

void
cs_parall_min(int             n,
              cs_datatype_t   datatype,
              void           *val);

#else

#define cs_parall_min(_n, _datatype, _val);

#endif

/*----------------------------------------------------------------------------
 * Broadcast values of a given datatype to all
 * default communicator processes.
 *
 * parameters:
 *   root_rank <-- rank from which to broadcast
 *   n         <-- number of values
 *   datatype  <-- matching Code_Saturne datatype
 *   val       <-- values to broadcast; input on root_rank,
 *                 output on others (size: n)
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Build a global array from each local array in each domain.
 *
 * Local arrays are appended in order of owning MPI rank.
 * The size of each local array may be different.
 *
 * Use of this function may be quite practical, but should be limited
 * to user functions, as it may limit scalability (especially as regards
 * memory usage).
 *
 * parameters:
 *   n_elts   <-- size of the local array
 *   n_g_elts <-- size of the global array
 *   array    <-- local array (size: n_elts)
 *   g_array  --> global array  (size: n_g_elts)
 *----------------------------------------------------------------------------*/

void
cs_parall_allgather_r(int        n_elts,
                      int        n_g_elts,
                      cs_real_t  array[],
                      cs_real_t  g_array[]);

/*----------------------------------------------------------------------------
 * Build an ordered global array from each local array in each domain.
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
 * parameters:
 *   n_elts   <-- number of local elements
 *   n_g_elts <-- number of global elements
 *   stride   <-- number of values per element
 *   o_key    <-- ordering key (coordinate) value per element
 *   array    <-- local array (size: n_elts*stride)
 *   g_array  --> global array  (size: n_g_elts*stride)
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Maximum value of a real and the value of related array on all
 * default communicator processes.
 *
 * parameters:
 *   n            <-- size of the related array
 *   max          <-> local max in, global max out
 *   max_loc_vals <-> array values at location of local max in,
 *                    and at location of global max out
 *----------------------------------------------------------------------------*/

void
cs_parall_max_loc_vals(int         n,
                       cs_real_t  *max,
                       cs_real_t   max_loc_vals[]);

/*----------------------------------------------------------------------------
 * Minimum value of a real and the value of related array on all
 * default communicator processes.
 *
 * parameters:
 *   n            <-- size of the related array
 *   min          <-> local min in, global min out
 *   min_loc_vals <-> array values at location of local min in,
 *                    and at location of global min out
 *----------------------------------------------------------------------------*/

void
cs_parall_min_loc_vals(int         n,
                       cs_real_t  *min,
                       cs_real_t   min_loc_vals[]);

/*----------------------------------------------------------------------------
 * Given an (id, rank, value) tuple, return the local id and rank
 * corresponding to the global minimum value.
 *
 * parameters:
 *   elt_id  <-> element id for which the value is the smallest
 *               (local in, global out)
 *   rank_id <-> rank id for which the value is the smallest
 *               (local in, global out)
 *   val     <-- associated local minimum value
 *----------------------------------------------------------------------------*/

void
cs_parall_min_id_rank_r(cs_lnum_t  *elt_id,
                        int        *rank_id,
                        cs_real_t   dis2mn);

/*----------------------------------------------------------------------------
 * Return minimum recommended scatter or gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * returns:
 *   minimum recommended part to block or gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

size_t
cs_parall_get_min_coll_buf_size(void);

/*----------------------------------------------------------------------------
 * Define minimum recommended gather buffer size.
 *
 * This is used by some internal part to block or scatter/gather algorithms,
 * so as to allow I/O buffer size tuning.
 *
 * parameters:
 *   minimum recommended part to block or gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

void
cs_parall_set_min_coll_buf_size(size_t buffer_size);

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
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARALL_H__ */
