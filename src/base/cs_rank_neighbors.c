/*============================================================================
 * Management of parallel rank neighbors.
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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_block_dist.h"
#include "cs_crystal_router.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rank_neighbors.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_rank_neighbors.c
        Management of parallel rank neighbors.

  Algorithm names are based upon \cite Hoefler:2010 and \cite Fox:1988 .
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/


/*=============================================================================
 * Local type definitions
 *============================================================================*/

#if defined(HAVE_MPI)


#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Exchange type names */

const char  *cs_rank_neighbors_exchange_name[]
= {"Personalized Exchange (PEX)",
   "Nonblocking Consensus (NBX)",
   "Crystal Router"};

/* Default synchronization mode */

static cs_rank_neighbors_exchange_t
_exchange_type = CS_RANK_NEIGHBORS_PEX;

#if defined(HAVE_MPI)

/* Call counters and timers (0: creation; 1: count/index; 2-4: communication
   based on method) */

static size_t              _rank_neighbors_calls[] = {0, 0, 0, 0, 0};
static cs_timer_counter_t  _rank_neighbors_timer[5];

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Descend binary tree for the sorting of an integer (ranks) array.
 *
 * parameters:
 *   level     <-- level of the binary tree to descend
 *   n         <-- number of entities in the binary tree to descend
 *   val       <-> integer array to order
 *----------------------------------------------------------------------------*/

inline static void
_heapsort_int_descend_tree(size_t         level,
                           const size_t   n,
                           int            val[])
{
  size_t  lv_cur;
  int v_save = val[level];

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n-1)
      if (val[lv_cur+1] > val[lv_cur]) lv_cur++;

    if (lv_cur >= n) break;

    if (v_save >= val[lv_cur]) break;

    val[level] = val[lv_cur];
    level = lv_cur;

  }

  val[level] = v_save;
}

/*----------------------------------------------------------------------------
 * Sort an array of local integers using heapsort
 *
 * parameters:
 *   val <-- array of values to sort
 *   n   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_heapsort_int(int           val[],
              const size_t  n)
{
  size_t i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2) ;
  do {
    i--;
    _heapsort_int_descend_tree(i, n, val);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n-1; i > 0 ; i--) {
    int v_save   = val[0];
    val[0] = val[i];
    val[i] = v_save;
    _heapsort_int_descend_tree(0, i, val);
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the sorting of an integer (ranks) array, with
 * an associated count.
 *
 * parameters:
 *   level     <-- level of the binary tree to descend
 *   n         <-- number of entities in the binary tree to descend
 *   val       <-> integer array to order
 *----------------------------------------------------------------------------*/

inline static void
_heapsort_int_descend_tree_count(size_t         level,
                                 const size_t   n,
                                 int            val[],
                                 cs_lnum_t      count[])
{
  size_t  lv_cur;
  int v_save = val[level];
  int c_save = count[level];

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n-1)
      if (val[lv_cur+1] > val[lv_cur]) lv_cur++;

    if (lv_cur >= n) break;

    if (v_save >= val[lv_cur]) break;

    val[level] = val[lv_cur];
    count[level] = count[lv_cur];
    level = lv_cur;

  }

  val[level] = v_save;
  count[level] = c_save;
}

/*----------------------------------------------------------------------------
 * Sort an array of local integers with associated count using heapsort
 *
 * parameters:
 *   val   <-- array of values to sort
 *   count <-- pre-allocated ordering table
 *   n     <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_heapsort_int_count(int           val[],
                    cs_lnum_t     count[],
                    const size_t  n)
{
  size_t i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2) ;
  do {
    i--;
    _heapsort_int_descend_tree_count(i, n, val, count);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n-1; i > 0 ; i--) {
    int v_save   = val[0];
    cs_lnum_t c_save   = count[0];
    val[0] = val[i];
    val[i] = v_save;
    count[0] = count[i];
    count[i] = c_save;
    _heapsort_int_descend_tree_count(0, i, val, count);
  }
}

/*----------------------------------------------------------------------------
 * Initialize rank neighbors using sort
 *
 * parameters:
 *   n        <-> rank neighbors structure
 *   n_elts   <-> number of elements
 *   elt_rank <-- array of element ranks
 *----------------------------------------------------------------------------*/

static void
_cs_rank_neighbors_init_sort(cs_rank_neighbors_t  *n,
                             size_t                n_elts,
                             const int             elt_rank[])
{
  BFT_MALLOC(n->rank, n_elts, int);

  /* Copy element ranks to sort tree, ignoring those whose values are
     are identical to the previous values (as some data  partitioning locality
     is usually quite good, this should lead to  a significantly smaller
     array to sort). */

  size_t n_elts_0 = 0;
  int r_prev = -1;
  for (size_t i = 0; i < n_elts; i++) {
    if (elt_rank[i] != r_prev) {
      n->rank[n_elts_0] = elt_rank[i];
      n_elts_0 += 1;
      r_prev = elt_rank[i];
    }
  }

  /* Now sort them */

  _heapsort_int(n->rank, n_elts_0);

  /* Remove duplicates */

  n->size = 0;
  r_prev = -1;
  for (size_t i = 0; i < n_elts_0; i++) {
    if (n->rank[i] != r_prev) {
      n->rank[n->size] = n->rank[i];
      n->size += 1;
      r_prev = n->rank[i];
    }
  }
  BFT_REALLOC(n->rank, n->size, int);
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a rank neighbors structure base on a list of element ranks
 *
 * \param[in]  n_elts    number of elements
 * \param[in]  elt_rank  element rank in
 *
 * \return  pointer to new rank neighborhood.
 */
/*----------------------------------------------------------------------------*/

cs_rank_neighbors_t *
cs_rank_neighbors_create(size_t     n_elts,
                         const int  elt_rank[])
{
  cs_rank_neighbors_t *n;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_rank_neighbors_calls[0] == 0)
    CS_TIMER_COUNTER_INIT(_rank_neighbors_timer[0]);

  /* Allocate structure */

  BFT_MALLOC(n, 1, cs_rank_neighbors_t);

  /* Create associated sub-structure */

  n->size = 0;
  n->rank = NULL;

  /* Determine neighboring ranks */

  _cs_rank_neighbors_init_sort(n, n_elts, elt_rank);

  /* timing */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_rank_neighbors_timer, &t0, &t1);
  _rank_neighbors_calls[0] += 1;

  return n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a rank neighborhood structure
 *
 * \param[in, out]  n   pointer to associated rank neighborhood
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_destroy(cs_rank_neighbors_t  **n)
{
  if (n != NULL) {

    cs_timer_t t0, t1;

    t0 = cs_timer_time();

    cs_rank_neighbors_t *_n = *n;

    BFT_FREE(_n->rank);
    BFT_FREE(*n);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(_rank_neighbors_timer, &t0, &t1);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Given a list of element ranks, determine the associated element
 *        rank indexes in a rank neighborhood structure.
 *
 * The elt_rank and elt_rank_index may be identical,
 * in which case it is updated.
 *
 * \param[in]   n               pointer to rank neighborhood structure
 * \param[in]   n_elts          number of elements
 * \param[in]   elt_rank        element rank (size: n_elts)
 * \param[out]  elt_rank_index  element rank index in neighborhood
 *                              (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_to_index(const cs_rank_neighbors_t  *n,
                           size_t                      n_elts,
                           const int                   elt_rank[],
                           int                        *elt_rank_index)
{
  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_rank_neighbors_calls[1] == 0)
    CS_TIMER_COUNTER_INIT(_rank_neighbors_timer[1]);

  /* Use binary search for each element of array */

  const size_t n_neighbors = n->size;
  const int  *rank = n->rank;

# pragma omp parallel if(n_elts > CS_THR_MIN)
  {
    int prev_rank = -1, prev_rank_index = -1;

#   pragma omp for
    for (size_t i = 0; i < n_elts; i++) {

      const int rank_1 = elt_rank[i];

      if (rank_1 == prev_rank)
        elt_rank_index[i] = prev_rank_index;

      else {

        /* use binary search */

        cs_lnum_t start_id = 0;
        cs_lnum_t end_id = n_neighbors - 1;
        cs_lnum_t mid_id = (end_id -start_id) / 2;

        while (start_id <= end_id) {
          int rank_cmp = rank[mid_id];
          if (rank_cmp < rank_1)
            start_id = mid_id + 1;
          else if (rank_cmp > rank_1)
            end_id = mid_id - 1;
          else
            break;
          mid_id = start_id + ((end_id -start_id) / 2);
        }

        if (rank[mid_id] == rank_1) {
          elt_rank_index[i] = mid_id;
          prev_rank_index =  mid_id;
        }
        else {
          elt_rank_index[i] = -1;
          prev_rank_index = -1;
        }

        prev_rank = rank_1;

      }

    }
  }

  /* Stop timer */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_rank_neighbors_timer + 1, &t0, &t1);
  _rank_neighbors_calls[1] += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Symmetrize a rank neighborhood structure.
 *
 * This is a collective operation ,which ensures that if rank i has j among
 * its neighbors, then j will also jave i among its neighbors.
 *
 * \param[in, out]  n     pointer to rank neighborhood structure
 * \param[in]       comm  associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_symmetrize(cs_rank_neighbors_t  *n,
                             MPI_Comm              comm)
{
  size_t      n_elts;
  cs_timer_t  t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_rank_neighbors_calls[2 + _exchange_type] == 0)
    CS_TIMER_COUNTER_INIT(_rank_neighbors_timer[2 + _exchange_type]);

  if (_exchange_type == CS_RANK_NEIGHBORS_PEX) {

    int comm_size;
    MPI_Comm_size(comm, &comm_size);

    int *sendbuf, *recvbuf;
    BFT_MALLOC(sendbuf, comm_size, int);
    BFT_MALLOC(recvbuf, comm_size, int);

    for (int i = 0; i < comm_size; i++)
      sendbuf[i] = 0;

    for (int i = 0; i < n->size; i++)
      sendbuf[n->rank[i]] = 1;

    MPI_Alltoall(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, comm);

    n_elts = 0;
    for (int i = 0; i < comm_size; i++) {
      if (recvbuf[i] != 0)
        n_elts++;
    }

    BFT_REALLOC(n->rank, n->size + n_elts, int);

    n_elts = 0;
    for (int i = 0; i < comm_size; i++) {
      if (recvbuf[i] != 0) {
        n->rank[n->size + n_elts] = i;
        n_elts++;
      }
    }

    n_elts += n->size;

    BFT_FREE(recvbuf);
    BFT_FREE(sendbuf);

  }

#if defined(HAVE_MPI_IBARRIER)

  /* Variant using nonblocking barrier */

  else if (_exchange_type == CS_RANK_NEIGHBORS_NBX) {

    size_t n_max_elts;

    int *sendbuf, *recvbuf;
    MPI_Request  *requests;

    BFT_MALLOC(sendbuf, n->size, int);
    BFT_MALLOC(requests, n->size, MPI_Request);

    n_elts = 0;
    n_max_elts = 16;
    BFT_MALLOC(recvbuf, n_max_elts, int);

    int tag = 0;

    for (int i = 0; i < n->size; i++) {
      sendbuf[i] = 1;
      MPI_Issend(sendbuf + i, 1, MPI_INT, n->rank[i], tag,
                 comm, requests + i);
    }

    MPI_Request  barrier_request;
    int flag;
    int barrier_done = 0, barrier_active = 0;

    while (!barrier_done) {

      MPI_Status status, status_recv;

      MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &flag, &status);
      if (flag) {
        if (n_elts >= n_max_elts) {
          n_max_elts *= 2;
          BFT_REALLOC(recvbuf, n_max_elts, int);
        }
        MPI_Recv(recvbuf + n_elts, 1, MPI_INT, status. MPI_SOURCE,
                 tag, comm, &status_recv);
        recvbuf[n_elts] = status.MPI_SOURCE;
        n_elts += 1;
      }

      if (!barrier_active) {

        MPI_Testall(n->size, requests, &flag, MPI_STATUSES_IGNORE);
        if (flag) {
          MPI_Ibarrier(comm, &barrier_request);
          barrier_active = 1;
        }

      }
      else
        MPI_Test(&barrier_request, &barrier_done, MPI_STATUS_IGNORE);

    }

    n_max_elts = n_elts;
    BFT_REALLOC(recvbuf, n_max_elts, int);

    BFT_REALLOC(n->rank, n->size + n_elts, int);
    for (size_t i = 0; i < n_elts; i++)
      n->rank[n->size + i] = recvbuf[i];

    n_elts += n->size;

    BFT_FREE(recvbuf);
    BFT_FREE(requests);
    BFT_FREE(sendbuf);

  }

#endif /* defined(HAVE_MPI_IBARRIER) */

  else { /* if (_exchange_type == CS_RANK_NEIGHBORS_CRYSTAL_ROUTER) */

    assert(_exchange_type == CS_RANK_NEIGHBORS_CRYSTAL_ROUTER);

    /* Create Crystal-router structure */

    int flags = CS_CRYSTAL_ROUTER_ADD_SRC_RANK;

    cs_crystal_router_t *cr = cs_crystal_router_create_s(n->size,
                                                         0,
                                                         CS_DATATYPE_NULL,
                                                         flags,
                                                         NULL,
                                                         NULL,
                                                         NULL,
                                                         n->rank,
                                                         comm);

    cs_crystal_router_exchange(cr);

    n_elts = cs_crystal_router_n_elts(cr);

    int *src_rank = NULL;
    cs_crystal_router_get_data(cr, &src_rank, NULL, NULL, NULL, NULL);

    BFT_REALLOC(n->rank, n->size + n_elts, int);
    for (size_t i = 0; i < n_elts; i++)
      n->rank[n->size + i] = src_rank[i];

    n_elts += n->size;

    BFT_FREE(src_rank);
    cs_crystal_router_destroy(&cr);

  }

  /* Sort by rank */

  _heapsort_int(n->rank, n_elts);

  /* Remove duplicates */

  n->size = 0;
  int r_prev = -1;
  for (size_t i = 0; i < n_elts; i++) {
    if (n->rank[i] != r_prev) {
      n->rank[n->size] = n->rank[i];
      n->size += 1;
      r_prev = n->rank[i];
    }
  }
  BFT_REALLOC(n->rank, n->size, int);

  /* Stop timer */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_rank_neighbors_timer + 2 + _exchange_type,
                            &t0, &t1);
  _rank_neighbors_calls[2 + _exchange_type] += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Given a list of element rank indexes, count occurences for
 *        a rank neighborhood structure.
 *
 * \param[in]   n               pointer to rank neighborhood structure
 * \param[in]   n_elts          number of elements
 * \param[in]   elt_rank_index  element rank index in neighborhood
 *                              (size: n_elts)
 * \param[out]  elt_rank_count  element rank count in neighborhood
 *                              (size: n->size)
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_count(const cs_rank_neighbors_t  *n,
                        size_t                      n_elts,
                        const int                  *elt_rank_index,
                        cs_lnum_t                  *elt_rank_count)
{
  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_rank_neighbors_calls[1] == 0)
    CS_TIMER_COUNTER_INIT(_rank_neighbors_timer[1]);

  for (int i = 0; i < n->size; i++)
    elt_rank_count[i] = 0;

  for (size_t i = 0; i < n_elts; i++)
    elt_rank_count[elt_rank_index[i]] += 1;

  /* Stop timer */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_rank_neighbors_timer + 1, &t0, &t1);
  _rank_neighbors_calls[1] += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange send and receive counts for rank neighborhoods.
 *
 * This allocates the n_recv ranks neighborhood structure and the
 * recv_count counts array, which the caller is responsible for freeing.
 *
 * \param[in]   n_send       pointer to rank neighborhood used for sending
 * \param[out]  n_recv       pointer to rank neighborhood used for receiving
 * \param[in]   send_count   pointer to rank neighborhood used for sending
 * \param[in]   recv_count   pointer to rank neighborhood used for sending
 * \param[in]   comm         associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_sync_count(const cs_rank_neighbors_t   *n_send,
                             cs_rank_neighbors_t        **n_recv,
                             const cs_lnum_t             *send_count,
                             cs_lnum_t                  **recv_count,
                             MPI_Comm                     comm)
{
  cs_rank_neighbors_sync_count_m(n_send,
                                 n_recv,
                                 send_count,
                                 recv_count,
                                 _exchange_type,
                                 comm);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange send and receive counts for rank neighborhoods,
 *        using a given method.
 *
 * This allocates the n_recv ranks neighborhood structure and the
 * recv_count counts array, which the caller is responsible for freeing.
 *
 * \param[in]   n_send         pointer to rank neighborhood used for sending
 * \param[out]  n_recv         pointer to rank neighborhood used for receiving
 * \param[in]   send_count     pointer to rank neighborhood used for sending
 * \param[in]   recv_count     pointer to rank neighborhood used for sending
 * \param[in]   exchange_type  exchange type
 * \param[in]   comm           associated communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_sync_count_m(const cs_rank_neighbors_t      *n_send,
                               cs_rank_neighbors_t           **n_recv,
                               const cs_lnum_t                *send_count,
                               cs_lnum_t                     **recv_count,
                               cs_rank_neighbors_exchange_t    exchange_type,
                               MPI_Comm                        comm)
{
  cs_timer_t  t0, t1;

  t0 = cs_timer_time();

  /* Allocate structure */

  cs_rank_neighbors_t  *_n_recv = NULL;
  cs_lnum_t            *_recv_count = NULL;

  BFT_MALLOC(_n_recv, 1, cs_rank_neighbors_t);
  _n_recv->rank = NULL;

  /* Fallback if nonblocking barrier is not available */
#if !defined(HAVE_MPI_IBARRIER)
  if (exchange_type == CS_RANK_NEIGHBORS_NBX)
    exchange_type == _exchange_type;
#endif

  /* Initialize timers if required */

  if (_rank_neighbors_calls[2 + _exchange_type] == 0)
    CS_TIMER_COUNTER_INIT(_rank_neighbors_timer[2 + exchange_type]);

  /* Exchange based on variant */

  switch(exchange_type) {

  case CS_RANK_NEIGHBORS_PEX:
    {

      int comm_size;
      MPI_Comm_size(comm, &comm_size);

      cs_lnum_t *sendbuf, *recvbuf;
      BFT_MALLOC(sendbuf, comm_size, cs_lnum_t);
      BFT_MALLOC(recvbuf, comm_size, cs_lnum_t);

      for (int i = 0; i < comm_size; i++)
        sendbuf[i] = 0;

      for (int i = 0; i < n_send->size; i++)
        sendbuf[n_send->rank[i]] = send_count[i];

      MPI_Alltoall(sendbuf, 1, CS_MPI_LNUM, recvbuf, 1, CS_MPI_LNUM, comm);

      _n_recv->size = 0;
      for (int i = 0; i < comm_size; i++) {
        if (recvbuf[i] != 0)
          _n_recv->size++;
      }

      BFT_MALLOC(_n_recv->rank, _n_recv->size, int);
      BFT_MALLOC(_recv_count, _n_recv->size, cs_lnum_t);

      size_t n_elts = 0;
      for (int i = 0; i < comm_size; i++) {
        if (recvbuf[i] != 0) {
          _n_recv->rank[n_elts] = i;
          _recv_count[n_elts] = recvbuf[i];
          n_elts++;
        }
      }

      BFT_FREE(recvbuf);
      BFT_FREE(sendbuf);

    }
    break;

#if defined(HAVE_MPI_IBARRIER)

  /* Variant using nonblocking barrier */

  case CS_RANK_NEIGHBORS_NBX:
    {

      size_t n_elts = 0;
      size_t n_max_elts = 16;

      MPI_Request  *requests;

      BFT_MALLOC(requests, n_send->size, MPI_Request);

      BFT_MALLOC(_n_recv->rank, n_max_elts, int);
      BFT_MALLOC(_recv_count, n_max_elts, cs_lnum_t);

      int tag = 0;

      for (int i = 0; i < n_send->size; i++) {
        MPI_Issend(send_count + i, 1, CS_MPI_LNUM, n_send->rank[i], tag,
                   comm, requests + i);
      }

      MPI_Request  barrier_request;
      int flag;
      int barrier_done = 0, barrier_active = 0;

      while (!barrier_done) {

        MPI_Status status, status_recv;

        MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &flag, &status);
        if (flag) {
          if (n_elts >= n_max_elts) {
            n_max_elts *= 2;
            BFT_REALLOC(_n_recv->rank, n_max_elts, int);
            BFT_REALLOC(_recv_count, n_max_elts, cs_lnum_t);
          }
          MPI_Recv(_recv_count + n_elts, 1, CS_MPI_LNUM, status.MPI_SOURCE,
                   tag, comm, &status_recv);
          _n_recv->rank[n_elts] = status.MPI_SOURCE;
          n_elts += 1;
        }

        if (!barrier_active) {

          MPI_Testall(n_send->size, requests, &flag, MPI_STATUSES_IGNORE);
          if (flag) {
            MPI_Ibarrier(comm, &barrier_request);
            barrier_active = 1;
          }

        }
        else
          MPI_Test(&barrier_request, &barrier_done, MPI_STATUS_IGNORE);

      }

      _n_recv->size = n_elts;
      BFT_REALLOC(_n_recv->rank, _n_recv->size, int);
      BFT_REALLOC(_recv_count, _n_recv->size, cs_lnum_t);

      _heapsort_int_count(_n_recv->rank, _recv_count, _n_recv->size);

      BFT_FREE(requests);

    }
    break;

#endif /* defined(HAVE_MPI_IBARRIER) */

  case CS_RANK_NEIGHBORS_CRYSTAL_ROUTER:
    {

      /* Create Crystal-router structure */

      int flags = CS_CRYSTAL_ROUTER_ADD_SRC_RANK;

      cs_crystal_router_t *cr = cs_crystal_router_create_s(n_send->size,
                                                           1,
                                                           CS_LNUM_TYPE,
                                                           flags,
                                                           send_count,
                                                           NULL,
                                                           NULL,
                                                           n_send->rank,
                                                           comm);

      cs_crystal_router_exchange(cr);

      _n_recv->size = cs_crystal_router_n_elts(cr);

      _n_recv->rank = NULL;

      void *recv_count_p = NULL;

      cs_crystal_router_get_data(cr,
                                 &(_n_recv->rank),
                                 NULL,
                                 NULL,
                                 NULL,
                                 &recv_count_p);

      _recv_count = recv_count_p;

      cs_crystal_router_destroy(&cr);

      _heapsort_int_count(_n_recv->rank, _recv_count, _n_recv->size);

    }
    break;

  default:
    assert(0);
    BFT_FREE(_n_recv);

  }

  /* Set return values */

  *n_recv     = _n_recv;
  *recv_count = _recv_count;

  /* Stop timer */

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_rank_neighbors_timer + 2 + _exchange_type,
                            &t0, &t1);
  _rank_neighbors_calls[2 + _exchange_type] += 1;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of rank neighbors collective algorithm choice.
 *
 * \return  current type of rank neighbors collective algorithm choice
 */
/*----------------------------------------------------------------------------*/

cs_rank_neighbors_exchange_t
cs_rank_neighbors_get_exchange_type(void)
{
  return _exchange_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of rank neighbors collective algorithm choice.
 *
 * \param  t  type of rank neighbors collective algorithm choice
 */
/*----------------------------------------------------------------------------*/

void
cs_rank_neighbors_set_exchange_type(cs_rank_neighbors_exchange_t  t)
{
#if !defined(HAVE_MPI_IBARRIER)
  if (t == CS_RANK_NEIGHBORS_NBX)
    bft_printf
      (_("Warning: The %s (%s)\n"
         "         exchange type is not available with the "
         "current MPI libary.\n"),
       "CS_RANK_NEIGHBORS_NBX", cs_rank_neighbors_exchange_name[t]);
  else
    _exchange_type = t;
#else
  _exchange_type = t;
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
