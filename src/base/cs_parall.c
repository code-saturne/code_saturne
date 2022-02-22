/*============================================================================
 * Functions dealing with parallelism
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "cs_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*! \fn inline static void cs_parall_counter(cs_gnum_t cpt[], const int n)
 *
 * \brief Sum values of a counter on all default communicator processes.
 *
 * \param[in, out]  cpt local counter in, global counter out (size: n)
 * \param[in]       n   number of values
 */

/*! \fn inline static void cs_parall_counter_max(cs_lnum_t cpt[], const int n)
 *
 * \brief Maximum values of a counter on all default communicator processes.
 *
 * \param[in, out]  cpt local counter in, global counter out (size: n)
 * \param[in]       n   number of values
 */

/*! \fn inline static void cs_parall_sum(int n, \
                                         cs_datatype_t datatype, \
                                         void *val)

 * \brief Sum values of a given datatype on all default communicator processes.
 *
 * \param[in]       n         number of values
 * \param[in]       datatype  matching code_saturne datatype
 * \param[in, out]  val       local sum in, global sum out (size: n)
 */

/*! \fn inline static void cs_parall_max(int n, \
                                         cs_datatype_t datatype, \
                                         void *val)
 *
 * \brief Maximum values of a given datatype on all
 *        default communicator processes.
 *
 * \param[in]       n         number of values
 * \param[in]       datatype  matching code_saturne datatype
 * \param[in, out]  val       local maximum in, global maximum out (size: n)
 */

/*! \fn inline static void cs_parall_min(int n, \
                                         cs_datatype_t datatype, \
                                         void *val)
 *
 *
 * \brief Minimum values of a given datatype on all
 *        default communicator processes.
 *
 * \param[in]       n         number of values
 * \param[in]       datatype  matching code_saturne datatype
 * \param[in, out]  val       local minimum in, global minimum out (size: n)
 */

/*! \fn inline static void cs_parall_bcast(int             root_rank, \
                                           int             n,         \
                                           cs_datatype_t   datatype,  \
                                           void           *val)
 *
 * \brief Broadcast values of a given datatype to all
 *        default communicator processes.
 *
 * \param[in]       root_rank  rank from which to broadcast
 * \param[in]       n          number of values
 * \param[in]       datatype   matching code_saturne datatype
 * \param[in, out]  val        values to broadcast; input on root_rank,
 *                             output on others (size: n)
 */

/*!
  \file cs_parall.c
        Utility functions dealing with parallelism.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

#define CS_PARALL_ARRAY_SIZE  500

typedef struct
{
  double  val;
  int     rank;
} _mpi_double_int_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

#if defined(HAVE_MPI)

/* Minimum recommended scatter/gather buffer size */

static size_t _cs_parall_min_coll_buf_size = 1024*1024;

#endif

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_parall_max_i(int  *max);

void
cs_f_parall_max_r(double  *max);

void
cs_f_parall_min_i(int  *min);

void
cs_f_parall_min_r(double  *min);

void
cs_f_parall_sum_i(int  *sum);

void
cs_f_parall_sum_r(double  *sum);

void
cs_f_parall_max_n_i(int   n,
                    int  *max);

void
cs_f_parall_max_n_r(int      n,
                    double  *max);

void
cs_f_parall_min_n_i(int   n,
                    int  *min);

void
cs_f_parall_min_n_r(int      n,
                    double  *min);

void
cs_f_parall_sum_n_i(int   n,
                    int  *sum);

void
cs_f_parall_sum_n_r(int      n,
                    double  *sum);

void
cs_f_parall_bcast_i(int   root_rank,
                    int  *val);

void
cs_f_parall_bcast_r(int      root_rank,
                    double  *val);

void
cs_f_parall_bcast_n_i(int   root_rank,
                      int   n,
                      int  *val);

void
cs_f_parall_bcast_n_r(int      root_rank,
                      int      n,
                      double  *val);

void
cs_f_parall_barrier(void);

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the array distribution over all ranks on a given root rank.
 *        The size of each local array may be different.
 *
 * \param[in]   root_rank  rank which stores the global array
 * \param[in]   n_elts     size of the local array
 * \param[in]   n_g_elts   size of the global array
 * \param[out]  count      number of elements in each rank
 * \param[out]  shift      shift to access each local array related to a rank
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

static void
_get_array_distribution(int     root_rank,
                        int     n_elts,
                        int     n_g_elts,
                        int    *p_count[],
                        int    *p_shift[])
{
  int  *count = NULL;
  int  *shift = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  assert(sizeof(double) == sizeof(cs_real_t));

  BFT_MALLOC(count, n_ranks, int);
  BFT_MALLOC(shift, n_ranks, int);

  MPI_Gather(&n_elts, 1, MPI_INT,
             count, 1, MPI_INT, root_rank, cs_glob_mpi_comm);

  shift[0] = 0;
  for (cs_lnum_t i = 1; i < n_ranks; i++)
    shift[i] = shift[i-1] + count[i-1];

  if (cs_glob_rank_id == root_rank)
    if (n_g_elts != (shift[n_ranks - 1] + count[n_ranks - 1]))
      bft_error(__FILE__, __LINE__, 0,
                _("Incorrect arguments to %s:\n"
                  "  sum of arg. 1 (n_elts) on ranks "
                  "is not equal to arg. 2 (n_g_elts)."), __func__);

  /* Return pointers */
  *p_count = count;
  *p_shift = shift;
}

#endif  /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Call MPI_Allreduce for a given code_saturne datatype and MPI
 * operation on all default communicator processes.
 *
 * parameters:
 *   n         <-- number of values
 *   datatype  <-- matching code_saturne datatype
 *   operation <-- MPI operation
 *   val       <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

static void
_cs_parall_allreduce(int             n,
                     cs_datatype_t   datatype,
                     MPI_Op          operation,
                     void           *val)
{
  int             data_size = n*cs_datatype_size[datatype];
  unsigned char  *locval;
  unsigned char  _locval[256];

  if (data_size > 256)
    BFT_MALLOC(locval, data_size, unsigned char);
  else
    locval = _locval;

  memcpy(locval, val, data_size);

  MPI_Allreduce(locval, val, n, cs_datatype_to_mpi[datatype], operation,
                cs_glob_mpi_comm);

  if (locval != _locval)
    BFT_FREE(locval);
}

#endif

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the global maximum of an integer in case of parallelism
 *
 * parameters:
 *   max <->  input = local max; output = global max
 *----------------------------------------------------------------------------*/

void
cs_f_parall_max_i(int  *max)
{
#if defined(HAVE_MPI)

  int  global_max;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (max, &global_max, 1, MPI_INT, MPI_MAX,
                 cs_glob_mpi_comm);

  *max = global_max;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global minimum of a real in case of parallelism
 *
 * parameters:
 *   min <->  input = local min; output = global min
 *----------------------------------------------------------------------------*/

void
cs_f_parall_min_r(double  *min)
{
#if defined(HAVE_MPI)

  double  global_min;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (min, &global_min, 1, MPI_DOUBLE, MPI_MIN,
                 cs_glob_mpi_comm);

  *min = global_min;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global minimum of an integer in case of parallelism
 *
 * parameters:
 *   min <->  input = local min; output = global min
 *----------------------------------------------------------------------------*/

void
cs_f_parall_min_i(int  *min)
{
#if defined(HAVE_MPI)

  int  global_min;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (min, &global_min, 1, MPI_INT, MPI_MIN,
                 cs_glob_mpi_comm);

  *min = global_min;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global maximum of a real in case of parallelism
 *
 * parameters:
 *   max <->  input = local max; output = global max
 *----------------------------------------------------------------------------*/

void
cs_f_parall_max_r(double  *max)
{
#if defined(HAVE_MPI)

  double  global_max;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                 cs_glob_mpi_comm);

  *max = global_max;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum of an integer in case of parallelism
 *
 * parameters:
 *   sum <->  input = local sum; output = global sum
 *----------------------------------------------------------------------------*/

void
cs_f_parall_sum_i(int  *sum)
{
#if defined(HAVE_MPI)

  int  global_sum;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (sum, &global_sum, 1, MPI_INT, MPI_SUM,
                 cs_glob_mpi_comm);

  *sum = global_sum;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global sum of a real in case of parallelism
 *
 * parameters:
 *   sum <->  input = local sum; output = global sum
 *----------------------------------------------------------------------------*/

void
cs_f_parall_sum_r(double  *sum)
{
#if defined(HAVE_MPI)

  double  global_sum;

  assert (sizeof (double) == sizeof (cs_real_t));

  MPI_Allreduce (sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                 cs_glob_mpi_comm);

  *sum = global_sum;

#endif
}

/*----------------------------------------------------------------------------
 * Compute the global maxima of an array of integers in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   max <-> input = local max; output = global max
 *----------------------------------------------------------------------------*/

void
cs_f_parall_max_n_i(int   n,
                    int  *max)
{
  cs_parall_max(n, CS_INT_TYPE, max);
}

/*----------------------------------------------------------------------------
 * Compute the global maxima of an array of reals in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   max <-> input = local max; output = global max
 *----------------------------------------------------------------------------*/

void
cs_f_parall_max_n_r(int      n,
                    double  *max)
{
  cs_parall_max(n, CS_DOUBLE, max);
}

/*----------------------------------------------------------------------------
 * Compute the global minima of an array of integers in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   min <-> input = local min; output = global min
 *----------------------------------------------------------------------------*/

void
cs_f_parall_min_n_i(int   n,
                    int  *min)
{
  cs_parall_min(n, CS_INT_TYPE, min);
}

/*----------------------------------------------------------------------------
 * Compute the global minima of an array of reals in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   min <-> input = local min; output = global min
 *----------------------------------------------------------------------------*/

void
cs_f_parall_min_n_r(int      n,
                    double  *min)
{
  cs_parall_min(n, CS_DOUBLE, min);
}

/*----------------------------------------------------------------------------
 * Compute the global sums of an array of integers in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   sum <-> input = local sum; output = global sum
 *----------------------------------------------------------------------------*/

void
cs_f_parall_sum_n_i(int   n,
                    int  *sum)
{
  cs_parall_sum(n, CS_INT_TYPE, sum);
}

/*----------------------------------------------------------------------------
 * Compute the global sums of an array of reals in case of parallelism
 *
 * parameters:
 *   n   <-- array size
 *   sum <-> input = local sum; output = global sum
 *----------------------------------------------------------------------------*/

void
cs_f_parall_sum_n_r(int      n,
                    double  *sum)
{
  cs_parall_sum(n, CS_DOUBLE, sum);
}

/*----------------------------------------------------------------------------
 * Broadcast an integer to all ranks
 *
 * parameters:
 *   root_rank <-- id of root rank
 *   val       <-> input = local value; output = global value
 *----------------------------------------------------------------------------*/

void
cs_f_parall_bcast_i(int   root_rank,
                    int  *val)
{
  cs_parall_bcast(root_rank, 1, CS_INT_TYPE, val);
}

/*----------------------------------------------------------------------------
 * Broadcast a double to all ranks
 *
 * parameters:
 *   root_rank <-- id of root rank
 *   val       <-> input = local value; output = global value
 *----------------------------------------------------------------------------*/

void
cs_f_parall_bcast_r(int      root_rank,
                    double  *val)
{
  cs_parall_bcast(root_rank, 1, CS_DOUBLE, val);
}

/*----------------------------------------------------------------------------
 * Broadcast an array of integers to all ranks
 *
 * parameters:
 *   root_rank <-- id of root rank
 *   n         <-- array size
 *   val       <-> input = local value; output = global value
 *----------------------------------------------------------------------------*/

void
cs_f_parall_bcast_n_i(int   root_rank,
                      int   n,
                      int  *val)
{
  cs_parall_bcast(root_rank, n, CS_INT_TYPE, val);
}

/*----------------------------------------------------------------------------
 * Broadcast an array of doubles to all ranks
 *
 * parameters:
 *   root_rank <-- id of root rank
 *   n         <-- array size
 *   val       <-> input = local value; output = global value
 *----------------------------------------------------------------------------*/

void
cs_f_parall_bcast_n_r(int      root_rank,
                      int      n,
                      double  *val)
{
  cs_parall_bcast(root_rank, n, CS_DOUBLE, val);
}

/*----------------------------------------------------------------------------
 * Call a barrier in case of parallelism
 *
 * This function should not be necessary in production code,
 * but it may be useful for debugging purposes.
 *----------------------------------------------------------------------------*/

void
cs_f_parall_barrier(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_counter(cs_gnum_t   cpt[],
                  const int   n)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, CS_GNUM_TYPE, MPI_SUM, cpt);
}

#endif

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_counter_max(cs_lnum_t   cpt[],
                      const int   n)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, CS_LNUM_TYPE, MPI_MAX, cpt);
}

#endif

/*----------------------------------------------------------------------------
 * Sum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching code_saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_sum(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_SUM, val);
}

#endif

/*----------------------------------------------------------------------------
 * Maximum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching code_saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_max(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_MAX, val);
}

#endif

/*----------------------------------------------------------------------------
 * Minimum values of a given datatype on all default communicator processes.
 *
 * parameters:
 *   n        <-- number of values
 *   datatype <-- matching code_saturne datatype
 *   val      <-> local value  input, global value output (array)
 *----------------------------------------------------------------------------*/

#if !defined(HAVE_MPI_IN_PLACE) && defined(HAVE_MPI)

void
cs_parall_min(int             n,
              cs_datatype_t   datatype,
              void           *val)
{
  if (cs_glob_n_ranks > 1)
    _cs_parall_allreduce(n, datatype, MPI_MIN, val);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Maximum value of a real and the value of related array on all
 * default communicator processes.
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
                       cs_real_t   max_loc_vals[])
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    _mpi_double_int_t  val_in, val_max;

    val_in.val  = *max;
    val_in.rank = cs_glob_rank_id;

    MPI_Allreduce(&val_in, &val_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                  cs_glob_mpi_comm);

    *max = val_max.val;

    MPI_Bcast(max_loc_vals, n, CS_MPI_REAL, val_max.rank, cs_glob_mpi_comm);

  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Minimum value of a real and the value of related array on all
 * default communicator processes.
 *
 * \param[in]       n             size of the related array
 * \param[in, out]  min           local min in, global min out
 * \param[in, out]  min_loc_vals  array values at location of local min in,
 *                                and at location of global min out
 */
/*----------------------------------------------------------------------------*/

void
cs_parall_min_loc_vals(int         n,
                       cs_real_t  *min,
                       cs_real_t   min_loc_vals[])
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    _mpi_double_int_t  val_in, val_min;

    val_in.val  = *min;
    val_in.rank = cs_glob_rank_id;

    MPI_Allreduce(&val_in, &val_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  cs_glob_mpi_comm);

    *min = val_min.val;

    MPI_Bcast(min_loc_vals, n, CS_MPI_REAL, val_min.rank, cs_glob_mpi_comm);

  }

#endif
}

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
                        cs_real_t   val)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_lnum_t buf[2];

    _mpi_double_int_t  val_in, val_min;

    assert(sizeof(double) == sizeof(cs_real_t));

    val_in.val  = val;
    val_in.rank = cs_glob_rank_id;

    MPI_Allreduce(&val_in, &val_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  cs_glob_mpi_comm);

    *rank_id = cs_glob_rank_id;

    buf[0] = *elt_id;
    buf[1] = *rank_id;

    MPI_Bcast(buf, 2, CS_MPI_LNUM, val_min.rank, cs_glob_mpi_comm);

    *elt_id = buf[0];
    *rank_id = buf[1];

  }

#endif
}

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
                      cs_real_t  g_array[])
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int  i;
    int  *count = NULL;
    int  *shift = NULL;

    const int  n_domains = cs_glob_n_ranks;

    assert(sizeof(double) == sizeof(cs_real_t));

    BFT_MALLOC(count, n_domains, int);
    BFT_MALLOC(shift, n_domains, int);

    MPI_Allgather(&n_elts, 1, MPI_INT, count, 1, MPI_INT,
                  cs_glob_mpi_comm);

    shift[0] = 0;
    for (i = 1; i < n_domains; i++)
      shift[i] = shift[i-1] + count[i-1];

    if (n_g_elts != (shift[n_domains - 1] + count[n_domains - 1]))
      bft_error(__FILE__, __LINE__, 0,
                _("Incorrect arguments to %s:\n"
                  "  sum of arg. 1 (n_elts) on ranks "
                  "is not equal to arg. 2 (n_g_elts)."),
                __func__);

    MPI_Allgatherv(array, n_elts, CS_MPI_REAL,
                   g_array, count, shift, CS_MPI_REAL, cs_glob_mpi_comm);

    BFT_FREE(count);
    BFT_FREE(shift);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    assert(n_elts == n_g_elts);

    for (int i = 0; i < n_elts; i++)
      g_array[i] = array[i];

  }
}

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
                              cs_real_t  g_array[])
{
  cs_lnum_t  *order;
  cs_real_t  *g_o_key;

  BFT_MALLOC(g_o_key, n_g_elts, cs_real_t);
  BFT_MALLOC(order, n_g_elts, cs_lnum_t);

  cs_parall_allgather_r(n_elts, n_g_elts, o_key, g_o_key);
  cs_parall_allgather_r(n_elts, n_g_elts*stride, array, g_array);

  cs_order_real_allocated(NULL, g_o_key, order, n_g_elts);
  cs_order_reorder_data(n_g_elts, sizeof(cs_real_t)*stride, order, g_array);

  BFT_FREE(order);
  BFT_FREE(g_o_key);
}

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
                   cs_real_t         g_array[])
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    /* Sanity check */

    if (cs_glob_rank_id == root_rank && g_array == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Global array is not allocated on the root_rank %d\n"),
                __func__, root_rank);

    int  *count = NULL, *shift = NULL;

    _get_array_distribution(root_rank, n_elts, n_g_elts, &count, &shift);

    MPI_Gatherv(array, n_elts, CS_MPI_REAL,
                g_array, count, shift, CS_MPI_REAL,
                root_rank, cs_glob_mpi_comm);

    BFT_FREE(count);
    BFT_FREE(shift);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    assert(n_elts == n_g_elts);
    assert(g_array != NULL && n_g_elts > 0);
    for (int i = 0; i < n_elts; i++)
      g_array[i] = array[i];

  }
}

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
                           cs_real_t  g_array[])
{
  cs_lnum_t  *order = NULL;
  cs_real_t  *g_o_key = NULL;

  if (cs_glob_rank_id == root_rank) {
    BFT_MALLOC(g_o_key, n_g_elts, cs_real_t);
    BFT_MALLOC(order, n_g_elts, cs_lnum_t);
  }

  cs_parall_gather_r(root_rank, n_elts, n_g_elts, o_key, g_o_key);
  cs_parall_gather_r(root_rank, n_elts, n_g_elts*stride, array, g_array);

  if (cs_glob_rank_id == root_rank) {
    cs_order_real_allocated(NULL, g_o_key, order, n_g_elts);
    cs_order_reorder_data(n_g_elts, sizeof(cs_real_t)*stride, order, g_array);
  }

  BFT_FREE(order);
  BFT_FREE(g_o_key);
}

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
                    cs_real_t         array[])
{
  assert(array != NULL && n_elts > 0);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    /* Sanity check */

    if (cs_glob_rank_id == root_rank && g_array == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Global array is not allocated on the root_rank %d\n"),
                __func__, root_rank);

    int  *count = NULL, *shift = NULL;

    _get_array_distribution(root_rank, n_elts, n_g_elts, &count, &shift);

    MPI_Scatterv(g_array, count, shift, CS_MPI_REAL,
                 array, n_elts, CS_MPI_REAL, root_rank, cs_glob_mpi_comm);

    BFT_FREE(count);
    BFT_FREE(shift);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    assert(n_elts == n_g_elts);
    assert(g_array != NULL && n_g_elts > 0);
    for (int i = 0; i < n_elts; i++)
      array[i] = g_array[i];

  }
}

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
                   float           g_array[])
{
  assert(array != NULL && n_elts > 0);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    /* Sanity check */

    if (cs_glob_rank_id == root_rank && g_array == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Global array is not allocated on the root_rank %d\n"),
                __func__, root_rank);

    int  *count = NULL, *shift = NULL;

    _get_array_distribution(root_rank, n_elts, n_g_elts, &count, &shift);

    MPI_Gatherv(array, n_elts, MPI_FLOAT,
                g_array, count, shift, MPI_FLOAT, root_rank, cs_glob_mpi_comm);

    BFT_FREE(count);
    BFT_FREE(shift);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    assert(n_elts == n_g_elts);
    assert(g_array != NULL && n_g_elts > 0);
    for (int i = 0; i < n_elts; i++)
      g_array[i] = array[i];

  }
}

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
                    float         array[])
{
  assert(array != NULL && n_elts > 0);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    /* Sanity check */

    if (cs_glob_rank_id == root_rank && g_array == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Global array is not allocated on the root_rank %d\n"),
                __func__, root_rank);

    int  *count = NULL, *shift = NULL;

    _get_array_distribution(root_rank, n_elts, n_g_elts, &count, &shift);

    MPI_Scatterv(g_array, count, shift, MPI_FLOAT,
                 array, n_elts, MPI_FLOAT, root_rank, cs_glob_mpi_comm);

    BFT_FREE(count);
    BFT_FREE(shift);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    assert(n_elts == n_g_elts);
    assert(g_array != NULL && n_g_elts > 0);
    for (int i = 0; i < n_elts; i++)
      array[i] = g_array[i];

  }
}

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
cs_parall_get_min_coll_buf_size(void)
{
#if defined(HAVE_MPI)
  return _cs_parall_min_coll_buf_size;
#else
  return 0;
#endif
}

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
cs_parall_set_min_coll_buf_size(size_t buffer_size)
{
#if defined(HAVE_MPI)
  _cs_parall_min_coll_buf_size = buffer_size;
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
