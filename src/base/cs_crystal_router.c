/*============================================================================
 * Crystal Router parallel exchange algorithm based operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_assert.h"
#include "cs_block_dist.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_crystal_router.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_crystal_router.c
        Crystal Router parallel exchange algorithm based operations.

  The Crystal Router is described in \cite Fox:1988 .

  \typedef cs_crystal_router_t
        Opaque Crystal Router distribution structure

  \par Using flags
  \parblock

  Flags are defined as a sum (bitwise or) of constants, which may include
  \ref CS_CRYSTAL_ROUTER_USE_DEST_ID, \ref CS_CRYSTAL_ROUTER_ADD_SRC_ID,
  and \ref CS_CRYSTAL_ROUTER_ADD_SRC_RANK.

  \endparblock

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Data elements in crystal router */

typedef enum {

  CS_CRYSTAL_ROUTER_DEST_ID,
  CS_CRYSTAL_ROUTER_SRC_ID,
  CS_CRYSTAL_ROUTER_DATA

} cs_crystal_router_data_t;

/*=============================================================================
 * Local type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/* Crystal router management structure */

struct  _cs_crystal_router_t { /* Crystal router information */

  cs_datatype_t   datatype;          /* associated datatype */
  int             flags;             /* ordering and metadata flags */

  size_t          stride;            /* stride if strided, 0 otherwise */

  size_t          dest_id_shift;     /* starting byte for destination id */
  size_t          src_id_shift;      /* starting byte for source id */
  size_t          elt_shift;         /* starting byte for element data */

  size_t          elt_size;          /* element size */
  size_t          comp_size;         /* composite metadata + element size if
                                        strided, 0 otherwise */
  size_t          n_elts[2];
  size_t          buffer_size[2];
  unsigned char  *buffer[2];

  MPI_Comm        comm;              /* associated MPI communicator */
  MPI_Datatype    comp_type;         /* Associated MPI datatype */
  int             rank_id;           /* local rank id in comm */
  int             n_ranks;           /* comm size */

};

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Call counter and timers: 0: total, 1: communication */

static size_t              _cr_calls = 0;
static cs_timer_counter_t  _cr_timers[2];

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * First stage of creation for a crystal router for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   flags            <-- ordering and metadata flags
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static cs_crystal_router_t *
_crystal_create_meta_s(size_t         n_elts,
                       int            stride,
                       cs_datatype_t  datatype,
                       int            flags,
                       MPI_Comm       comm)
{
  int rank_id, n_ranks;
  cs_crystal_router_t *cr = NULL;

  size_t comp_size = 0;
  size_t elt_size = cs_datatype_size[datatype]*stride;
  size_t align_size = sizeof(cs_lnum_t);

  /* Ensure alignement on integer size at least */

  if (elt_size % sizeof(int) > 0)
    comp_size += sizeof(int) - elt_size % sizeof(int);

  /* Communicator info */

  MPI_Comm_rank(comm, &rank_id);
  MPI_Comm_size(comm, &n_ranks);

  /* Allocate structure */

  BFT_MALLOC(cr, 1, cs_crystal_router_t);

  cr->datatype = (stride > 0) ? datatype : CS_DATATYPE_NULL;

  cr->flags = flags;

  cr->stride = (stride > 0) ? stride : 1;
  cr->elt_size = elt_size;
  cr->comp_size = comp_size;
  cr->n_elts[0] = n_elts;
  cr->n_elts[1] = 0;

  cr->comm = comm;
  cr->rank_id = rank_id;
  cr->n_ranks = n_ranks;

  /* Compute data size and alignment */

  cr->dest_id_shift = sizeof(int);

  if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_RANK)
    cr->dest_id_shift += sizeof(int);

  if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID)
    cr->src_id_shift = cr->dest_id_shift + sizeof(cs_lnum_t);
  else
    cr->src_id_shift = cr->dest_id_shift;

  if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_ID)
    cr->elt_shift = cr->src_id_shift + cs_datatype_size[CS_LNUM_TYPE];
  else
    cr->elt_shift = cr->src_id_shift;

  if (cr->elt_shift % align_size)
    cr->elt_shift += align_size - (cr->elt_shift % align_size);

  cr->comp_size = cr->elt_shift + elt_size;

  if (elt_size % align_size)
    cr->comp_size += align_size - (elt_size % align_size);;

  /* Create associated MPI datatype */

  MPI_Type_contiguous(cr->comp_size, MPI_BYTE, &(cr->comp_type));
  MPI_Type_commit(&(cr->comp_type));

  /* Allocate buffers */

  cr->buffer_size[0] = n_elts*cr->comp_size;
  cr->buffer_size[1] = 0;
  BFT_MALLOC(cr->buffer[0], cr->buffer_size[0], unsigned char);
  memset(cr->buffer[0], 0, cr->buffer_size[0]);
  cr->buffer[1] = NULL;

  return cr;
}

/*----------------------------------------------------------------------------
 * Partition strided data for exchange with a crystal router.
 *
 * parameters:
 *   cr      <-> associated crystal router structure
 *   send_id <-> id of "send" buffer" (1 if low = buf0, 0 if low = buf1)
 *   cutoff  <-- cutoff rank
 *---------------------------------------------------------------------------*/

static void
_crystal_partition_strided(cs_crystal_router_t  *cr,
                           int                   send_id,
                           int                   cutoff)
{
  cs_lnum_t i;

  cs_lnum_t n0 = 0, n1 = 0;

  const cs_lnum_t n = cr->n_elts[0];
  const int id0 = (send_id + 1) % 2;
  const int id1 = send_id;
  const size_t comp_size = cr->comp_size;

  assert(send_id == 0 || send_id == 1);

  if (cr->buffer_size[1] < cr->buffer_size[0]) {
    cr->buffer_size[1] = cr->buffer_size[0];
    BFT_REALLOC(cr->buffer[1], cr->buffer_size[1], unsigned char);
  }

  if (id0 == 0) {
    for (i = 0; i < n; i++) {
      unsigned char *src = cr->buffer[0] + i*comp_size;
      int *r = (int *)src;
      if (r[0] >= cutoff)
        break;
    }
    n0 = i;
  }
  else {
    for (i = 0; i < n; i++) {
      unsigned char *src = cr->buffer[0] + i*comp_size;
      int *r = (int *)src;
      if (r[0] < cutoff)
        break;
    }
    n1 = i;
  }

  while (i < n) {
    unsigned char *src = cr->buffer[0] + i*comp_size;
    int *r = (int *)src;
    if (r[0] < cutoff) {
      unsigned char *dest = cr->buffer[id0] + n0*comp_size;
      memcpy(dest, src, comp_size);
      n0++;
    }
    else {
      unsigned char *dest = cr->buffer[id1] + n1*comp_size;
      memcpy(dest, src, comp_size);
      n1++;
    }
    i++;
  }

  assert(n0 + n1 == n);

  cr->n_elts[id0] = n0;
  cr->n_elts[id1] = n1;
}

/*----------------------------------------------------------------------------
 * Send and receive strided data with a crystal router for one stage.
 *
 * parameters:
 *   cr        <-> associated crystal router structure
 *   target    <-- base target rank to send/receive from
 *   n_recv    <-- number of ranks to receive from
 *---------------------------------------------------------------------------*/

static  void
_crystal_sendrecv_strided(cs_crystal_router_t  *cr,
                          int                   target,
                          int                   n_recv)
{
  cs_timer_t t0 = cs_timer_time();

  int send_size;
  cs_lnum_t i;
  size_t loc_size;
  uint64_t test_size;

  MPI_Status status[3];
  MPI_Request request[3] = {MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL};
  int recv_size[2] = {0, 0};

  /* Send message to target process */

  send_size = cr->n_elts[1];

  test_size = (uint64_t)cr->n_elts[1];
  if ((uint64_t)send_size != test_size)
    bft_error(__FILE__, __LINE__, 0,
              "Crystal router:"
              "  Message to send would have size too large for C int: %llu",
              (unsigned long long)test_size);

  MPI_Isend(&send_size, 1, MPI_INT, target, cr->rank_id,
            cr->comm, &request[0]);

  for (i = 0; i < n_recv; i++)
    MPI_Irecv(recv_size+i, 1, MPI_INT, target+i, target+i,
              cr->comm, request+i+1);

  MPI_Waitall(n_recv + 1, request, status);

  loc_size = cr->n_elts[0]*cr->comp_size;
  for (i=0; i < n_recv; i++)
    loc_size += recv_size[i]*cr->comp_size;

  if (loc_size > cr->buffer_size[0]) {
    cr->buffer_size[0] = loc_size;
    BFT_REALLOC(cr->buffer[0], cr->buffer_size[0], unsigned char);
  }

  MPI_Isend(cr->buffer[1], cr->n_elts[1], cr->comp_type,
            target, cr->rank_id, cr->comm, request);

  cr->n_elts[1] = 0;

  if (n_recv) {

    unsigned char *r_ptr = cr->buffer[0] + (cr->n_elts[0]*cr->comp_size);

    MPI_Irecv(r_ptr, recv_size[0], cr->comp_type,
              target, target, cr->comm, request+1);

    cr->n_elts[0] += recv_size[0];

    if (n_recv == 2) {

      r_ptr += recv_size[0]*cr->comp_size;

      MPI_Irecv(r_ptr, recv_size[1], cr->comp_type,
                target+1, target+1, cr->comm, request+2);

      cr->n_elts[0] += recv_size[1];

    }

  }

  MPI_Waitall(n_recv+1, request, status);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers + 1, &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Exchange strided data with a crystal router.
 *
 * Order of data from a same source rank is preserved.
 *
 * parameters:
 *   cr <-> associated crystal router structure
 *---------------------------------------------------------------------------*/

static  void
_crystal_exchange_strided(cs_crystal_router_t  *cr)
{
  int target = -1;
  int send_part = 0;
  int b_low = 0;
  int n_sub_ranks = cr->n_ranks;

  while (n_sub_ranks > 1) {

    int n_recv = 1;
    int n_low = n_sub_ranks / 2;
    int b_high = b_low + n_low;

    if (cr->rank_id < b_high) {
      target = cr->rank_id + n_low;
      if ((n_sub_ranks & 1) && (cr->rank_id == b_high - 1))
        n_recv = 2;
      send_part = 1;
    }
    else {
      target = cr->rank_id - n_low;
      if (target == b_high) {
        target--;
        n_recv = 0;
      }
      send_part = 0;
    }

    /* Partition data */

    _crystal_partition_strided(cr, send_part, b_high);

    /* Send message to target process */

    _crystal_sendrecv_strided(cr, target, n_recv);

    /* Ready for next exchange */

    if (cr->rank_id < b_high)
      n_sub_ranks = n_low;
    else {
      n_sub_ranks -= n_low;
      b_low = b_high;
    }

  }

  cr->n_elts[1] = 0;
  cr->buffer_size[1] = 0;
  BFT_FREE(cr->buffer[1]);
}

/*----------------------------------------------------------------------------
 * Compare strided elements from crystal router (qsort function).
 *
 * parameters:
 *   x <-> pointer to first element
 *   y <-> pointer to second element
 *
 * returns:
 *   -1 if x < y, 0 if x = y, or 1 if x > y
 *----------------------------------------------------------------------------*/

static int _compare_strided(const void *x, const void *y)
{
  int retval = 1;

  const int *c0 = x;
  const int *c1 = y;

  if (c0[1] < c1[1])
    retval = -1;

  /* For same rank, preserve relative positions, so compare pointer */

  else if (c0[1] == c1[1]) {
    if (c0 < c1)
      retval = -1;
    else if (c0 == c1)
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Sort strided crystal router data by source rank.
 *
 * parameters:
 *   cr      <-> associated crystal router structure
 *   send_id <-> id of "send" buffer" (1 if low = buf0, 0 if low = buf1)
 *   cutoff  <-- cutoff rank
 *---------------------------------------------------------------------------*/

static void
_crystal_sort_by_source_rank_strided(cs_crystal_router_t  *cr)
{
  const cs_lnum_t n = cr->n_elts[0];

  if (n > 0)
    qsort(cr->buffer[0], n, cr->comp_size, &_compare_strided);
}

/*----------------------------------------------------------------------------
 * Get data elements associated with a strided Crystal Router.
 *
 * This variant assumes no destination id is provided, so data is
 * copied in Crystal Router buffer order.
 *
 * parameters:
 *   <-- cr         pointer to associated Crystal Router
 *   --> src_rank   pointer to source rank array, or NULL
 *   --> src_id     pointer to source id array, or NULL
 *   --> data  pointer to (pointer to) destination data, or NULL
 *----------------------------------------------------------------------------*/

static void
_get_data_s(cs_crystal_router_t   *cr,
            int                   *src_rank,
            cs_lnum_t             *src_id,
            void                  *data)
{
  const size_t n_elts = cr->n_elts[0];

  /* Now extract data; work with chunks to amortize some tests
     while not multiplying variants */

  const size_t chunk_size = 64;
  const size_t n_chunks
    = (n_elts % chunk_size) ? n_elts/chunk_size + 1 : n_elts/chunk_size;

  unsigned char *cr_src_rank_p = cr->buffer[0] + sizeof(int);
  unsigned char *cr_src_id_p = cr->buffer[0] + cr->src_id_shift;
  unsigned char *cr_data_p = cr->buffer[0] + cr->elt_shift;

  unsigned char *src_rank_p = (unsigned char *)src_rank;
  unsigned char *src_id_p = (unsigned char *)src_id;
  unsigned char *data_p = (unsigned char *)data;

  for (size_t c_id = 0; c_id < n_chunks; c_id++) {

    size_t i0 = c_id*chunk_size;
    size_t i1 = (c_id+1)*chunk_size;
    if (i1 > n_elts)
      i1 = n_elts;

    /* Optional source rank */

    if (src_rank != NULL) {
      for (size_t i = i0; i < i1; i++) {
        for (size_t j = 0; j < sizeof(int); j++)
          src_rank_p[i*sizeof(int) + j]
            = cr_src_rank_p[i*cr->comp_size + j];
      }
    }

    /* Optional source id */

    if (src_id != NULL) {
      for (size_t i = i0; i < i1; i++) {
        for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
          src_id_p[i*sizeof(cs_lnum_t) + j]
            = cr_src_id_p[i*cr->comp_size + j];
      }
    }

    /* data */

    if (data != NULL) {
      for (size_t i = i0; i < i1; i++) {
        for (size_t j = 0; j < cr->elt_size; j++)
          data_p[i*cr->elt_size + j]
            = cr_data_p[i*cr->comp_size + j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Get data elements associated with a strided Crystal Router.
 *
 * This variant assumes a destination id provided in the Crystal Router,
 * which will be applied automatically to all extracted arrays, except the
 * dest_id array itself, which is always in buffer order.
 *
 * parameters:
 *   <-- cr         pointer to associated Crystal Router
 *   --> src_rank   pointer to source rank array, or NULL
 *   --> dest_id    pointer to destination id array, or NULL
 *   --> src_id     pointer to source id array, or NULL
 *   --> data  pointer to (pointer to) destination data, or NULL
 *----------------------------------------------------------------------------*/

static void
_get_data_s_with_dest_id(cs_crystal_router_t   *cr,
                         int                   *src_rank,
                         cs_lnum_t             *dest_id,
                         cs_lnum_t             *src_id,
                         void                  *data)
{
  cs_lnum_t *_dest_id = dest_id;

  const size_t n_elts = cr->n_elts[0];

  /* Now extract data; work with chunks to amortize some tests
     while not multiplying variants */

  const size_t chunk_size = 64;
  const size_t n_chunks
    = (n_elts % chunk_size) ? n_elts/chunk_size + 1 : n_elts/chunk_size;

  unsigned char *cr_src_rank_p = cr->buffer[0] + sizeof(int);
  unsigned char *cr_dest_id_p = cr->buffer[0] + cr->dest_id_shift;
  unsigned char *cr_src_id_p = cr->buffer[0] + cr->src_id_shift;
  unsigned char *cr_data_p = cr->buffer[0] + cr->elt_shift;

  /* If we do not have a dest_id buffer, use a local one */

  if (_dest_id == NULL) {
    BFT_MALLOC(_dest_id, n_elts, cs_lnum_t);
    cs_assert(cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID);
  }

  unsigned char *dest_id_p = (unsigned char *)_dest_id;
  unsigned char *src_rank_p = (unsigned char *)src_rank;
  unsigned char *src_id_p = (unsigned char *)src_id;
  unsigned char *data_p = (unsigned char *)data;

  for (size_t c_id = 0; c_id < n_chunks; c_id++) {

    size_t i0 = c_id*chunk_size;
    size_t i1 = (c_id+1)*chunk_size;
    if (i1 > n_elts)
      i1 = n_elts;

    /* Extract _dest_id first, to use if in other loops */

    if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
      for (size_t i = i0; i < i1; i++) {
        for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
          dest_id_p[i*sizeof(cs_lnum_t) + j]
            = cr_dest_id_p[i*cr->comp_size + j];
      }
    }

    /* Optional source rank */

    if (src_rank != NULL) {
      for (size_t i = i0; i < i1; i++) {
        cs_lnum_t e_id = _dest_id[i];
        for (size_t j = 0; j < sizeof(int); j++)
          src_rank_p[e_id*sizeof(int) + j]
            = cr_src_rank_p[i*cr->comp_size + j];
      }
    }

    /* Optional source id */

    if (src_id != NULL) {
      for (size_t i = i0; i < i1; i++) {
        cs_lnum_t e_id = _dest_id[i];
        for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
          src_id_p[e_id*sizeof(cs_lnum_t) + j]
            = cr_src_id_p[i*cr->comp_size + j];
      }
    }

    /* data */

    if (data != NULL) {
      for (size_t i = i0; i < i1; i++) {
        cs_lnum_t e_id = _dest_id[i];
        for (size_t j = 0; j < cr->elt_size; j++)
          data_p[e_id*cr->elt_size + j]
            = cr_data_p[i*cr->comp_size + j];
      }
    }

  }

  if (dest_id == NULL)
    BFT_FREE(_dest_id);
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Crystal Router for strided data.
 *
 * If the flags constant contains \ref CS_CRYSTAL_ROUTER_USE_DEST_ID,
 * data exchanged will be ordered by the array passed to the
 * \c dest_id argument. For \c n total values received on a rank
 * (as given by \ref cs_crystal_router_n_elts), those destination ids
 * must be in the range [0, \c n[.
 *
 * If the flags bit mask matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_ID,
 * source ids are added. If it matches \ref CS_CRYSTAL_ROUTER_ADD_SRC_RANK,
 * source rank metadata is added.
 *
 * \param[in]  n_elts            number of elements
 * \param[in]  stride            number of values per entity (interlaced)
 * \param[in]  datatype          type of data considered
 * \param[in]  flags             add destination id ?
 * \param[in]  elt               element values
 * \param[in]  dest_id           element destination id, or NULL
 * \param[in]  dest_rank         destination rank for each element
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new Crystal Router structure.
 */
/*----------------------------------------------------------------------------*/

cs_crystal_router_t *
cs_crystal_router_create_s(size_t            n_elts,
                           int               stride,
                           cs_datatype_t     datatype,
                           int               flags,
                           const void       *elt,
                           const cs_lnum_t  *dest_id,
                           const int         dest_rank[],
                           MPI_Comm          comm)
{
  cs_timer_t t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_cr_calls == 0) {
    for (int i = 0; i < 2; i++)
      CS_TIMER_COUNTER_INIT(_cr_timers[i]);
  }
  _cr_calls += 1;

  unsigned const char *_elt = elt;

  /* Allocate structure */

  cs_crystal_router_t *cr = _crystal_create_meta_s(n_elts,
                                                   stride,
                                                   datatype,
                                                   flags,
                                                   comm);
  /* Copy data and some metadata */

  const bool add_src_rank
    = (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_RANK) ? true : false;
  const bool add_dest_id
    = (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) ? true : false;
  const bool add_src_id
    = (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_ID) ? true : false;

  if (add_dest_id)
    cs_assert(dest_id != NULL || n_elts == 0);

  for (size_t i = 0; i < n_elts; i++) {

    int *pr = (int *)(cr->buffer[0] + i*cr->comp_size);
    unsigned char *pe = cr->buffer[0] + i*cr->comp_size + cr->elt_shift;
    unsigned const char *_psrc = _elt + i*cr->elt_size;

    pr[0] = dest_rank[i];

    if (add_src_rank)
      pr[1] = cr->rank_id;

    if (add_dest_id) {
      unsigned char *pi = cr->buffer[0] + i*cr->comp_size + cr->dest_id_shift;
      unsigned const char *_p_dest_id = (unsigned const char *)(dest_id + i);
      for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
        pi[j] = _p_dest_id[j];
    }

    if (add_src_id) {
      cs_lnum_t src_id = i;
      unsigned char *_src_id = (unsigned char *)(&src_id);
      unsigned char *pi = cr->buffer[0] + i*cr->comp_size + cr->src_id_shift;
      for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
        pi[j] = _src_id[j];
    }

    for (size_t j = 0; j < cr->elt_size; j++)
      pe[j] = _psrc[j];
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);

  return cr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a Crystal Router.
 *
 * \param[in, out]  cr   pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_destroy(cs_crystal_router_t  **cr)
{
  if (cr != NULL) {

    cs_timer_t t0 = cs_timer_time();

    if (cr != NULL) {
      cs_crystal_router_t *_cr = *cr;
      if (_cr->comp_type != MPI_BYTE)
        MPI_Type_free(&(_cr->comp_type));
      BFT_FREE(_cr->buffer[1]);
      BFT_FREE(_cr->buffer[0]);
      BFT_FREE(*cr);
    }

    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(_cr_timers, &t0, &t1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange data with a Crystal Router.
 *
 * Order of data from a same source rank is preserved.
 *
 * \param[in, out]  cr   pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_exchange(cs_crystal_router_t  *cr)
{
  cs_assert(cr != NULL);

  cs_timer_t t0 = cs_timer_time();

  /* Only strided version implemented so far */

  _crystal_exchange_strided(cr);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sort stride crystal router data by source rank.
 *
 * \param[in, out]  cr  pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_sort_by_source_rank(cs_crystal_router_t  *cr)
{
  cs_assert(cr != NULL);

  cs_timer_t t0 = cs_timer_time();

  /* Only strided version implemented so far */

  _crystal_sort_by_source_rank_strided(cr);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements associated with Crystal Router.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * \param[in]  cr  pointer to associated Crystal Router
 *
 * \return  number of elements associated with distributor.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_crystal_router_n_elts(const cs_crystal_router_t  *cr)
{
  cs_lnum_t retval = 0;

  if (cr != NULL)
    retval = cr->n_elts[0];

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get data elements associated with a Crystal Router.
 *
 * Depending on the creation options, some elements may be available or not.
 *
 * For each element type, a pointer to a pointer to a buffer is provided.
 * If the pointer to the buffer is non-null, the buffer must be of sufficient
 * size for the number of elements returned by \ref cs_crystal_router_n_elts.
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * Note also that if the destination id is provided in the Crystal Router,
 * it will be applied automatically to all extracted arrays, except the
 * dest_id array itself, which is always in receive order. If the Crystal
 * Router does not contain destination id info but the \c dest_id
 * argument points to a non-NULL value, the provided id will be used to
 * order extracted data. This allows saving the destination id on the receive
 * side, and not re-sending it (saving bandwidth) for subsequent calls
 * with a similar Crystal Router.
 *
 * \param[in]   cr          pointer to associated Crystal Router
 * \param[out]  src_rank    pointer to (pointer to) source rank array, or NULL
 * \param[out]  dest_id     pointer to (pointer to) destination id array,
 *                          or NULL
 * \param[out]  src_id      pointer to (pointer to) source id array, or NULL
 * \param[out]  dest_index  pointer to (pointer to) destination index, or NULL
 * \param[out]  dest_data   pointer to (pointer to) destination data, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_get_data(cs_crystal_router_t   *cr,
                           int                  **src_rank,
                           cs_lnum_t            **dest_id,
                           cs_lnum_t            **src_id,
                           cs_lnum_t            **dest_index,
                           void                 **dest_data)
{
  cs_timer_t t0 = cs_timer_time();

  int *_src_rank = NULL;
  cs_lnum_t *_dest_id = NULL;
  cs_lnum_t *_src_id = NULL;
  cs_lnum_t *_dest_index = NULL;
  unsigned char *_dest_data = NULL;

  size_t n_elts = cr->n_elts[0];

  /* Map and allocate arrays if needed */

  if (src_rank != NULL && cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_RANK) {
    _src_rank = *src_rank;
    if (_src_rank == NULL) {
      BFT_MALLOC(_src_rank, n_elts, int);
      *src_rank = _src_rank;
    }
  }

  if (dest_id != NULL) { /* May be extracted or provided from previous call */
    _dest_id = *dest_id;
    if (_dest_id == NULL && cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
      BFT_MALLOC(_dest_id, n_elts, cs_lnum_t);
      *dest_id = _dest_id;
    }
  }

  if (src_id != NULL && cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_ID) {
    _src_id = *src_id;
    if (_src_id == NULL) {
      BFT_MALLOC(_src_id, n_elts, cs_lnum_t);
      *src_id = _src_id;
    }
  }

  if (dest_index != NULL && cr->stride == 0) {
    _dest_index = *dest_index;
    if (_dest_index == NULL) {
      BFT_MALLOC(_dest_index, n_elts + 1, cs_lnum_t);
      _dest_index[0] = 0;
      *dest_index = _dest_index;
    }
  }

  if (dest_data != NULL && cr->stride > 0) {
    _dest_data = *dest_data;
    if (_dest_data == NULL) {
      BFT_MALLOC(_dest_data, n_elts*cr->elt_size, unsigned char);
      *dest_data = _dest_data;
    }
  }

  /* Now extract data */

  if (cr->stride > 0) {

    if (_dest_id != NULL || cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID)
      _get_data_s_with_dest_id(cr,
                               _src_rank,
                               _dest_id,
                               _src_id,
                               _dest_data);

    else
      _get_data_s(cr, _src_rank, _src_id, _dest_data);

  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log performance information relative to Crystal Router exchange.
 *
 * Call count is based on Crystal router creation calls, as the number
 * of exchange and data access calls should be identical.
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_log_finalize(void)
{
  if (_cr_calls <= 0 || cs_glob_n_ranks < 2)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nCrystal router: %llu %s\n"),
                (unsigned long long) _cr_calls, _("calls"));

#if defined(HAVE_MPI)
  double wtimes[2] = {_cr_timers[0].wall_nsec*1e-9,
                      _cr_timers[1].wall_nsec*1e-9};

  double mntimes[2], mxtimes[2], stimes[2];

  MPI_Reduce(wtimes, mntimes, 2, MPI_DOUBLE, MPI_MIN,
             0, cs_glob_mpi_comm);
  MPI_Reduce(wtimes, mxtimes, 2, MPI_DOUBLE, MPI_MAX,
             0, cs_glob_mpi_comm);
  MPI_Reduce(wtimes, stimes, 2, MPI_DOUBLE, MPI_SUM,
             0, cs_glob_mpi_comm);
  if (cs_glob_rank_id == 0) {
    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("                      mean           minimum        maximum\n"
         "  wall clock:    %12.5f s %12.5f s %12.5f s\n"
         "  communication: %12.5f s %12.5f s %12.5f s\n"),
       stimes[0]/cs_glob_n_ranks, mntimes[0], mxtimes[0],
       stimes[1]/cs_glob_n_ranks, mntimes[1], mxtimes[1]);
  }
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
