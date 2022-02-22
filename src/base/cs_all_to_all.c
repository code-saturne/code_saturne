/*============================================================================
 * All-to-all parallel data exchange.
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
#include "bft_mem_usage.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_assert.h"
#include "cs_block_dist.h"
#include "cs_crystal_router.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_rank_neighbors.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_all_to_all.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_all_to_all.c
        All-to-all parallel data exchange.

  \typedef cs_all_to_all_t
        Opaque all-to-all distribution structure

  \enum cs_all_to_all_type_t

  \brief All-to-all algorithm selection

  \var CS_ALL_TO_ALL_MPI_DEFAULT
       Use MPI_Alltoall and MPI_Alltoallv sequences

  \var CS_ALL_TO_ALL_CRYSTAL_ROUTER
       Use crystal router algorithm

  \paragraph all_to_all_flags Using flags
  \parblock

  Flags are defined as a sum (bitwise or) of constants, which may include
  \ref CS_ALL_TO_ALL_USE_DEST_ID, \ref CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
  \ref CS_ALL_TO_ALL_NO_REVERSE, and \ref CS_ALL_TO_ALL_NEED_SRC_RANK.

  \endparblock
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Note for future optimization:
   for cases where data to be sent is sorted by increasing destination rank,
   using communication schemes requiring ordering by ranks (such as
   CS_ALL_TO_ALL_MPI_DEFAULT, we should not need to build an extra array,
   but only to send the correct parts of the indexed list to the correct
   ranks, so cs_all_to_all_copy_indexed may internally use an extra array;
   this could be avoided by adding an additional send option or flag,
   or adding cs_all_to_copy_... variants for cases ordered by rank */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

typedef enum {

  CS_ALL_TO_ALL_TIME_TOTAL,
  CS_ALL_TO_ALL_TIME_METADATA,
  CS_ALL_TO_ALL_TIME_EXCHANGE

} cs_all_to_all_timer_t;

/* Base structure for MPI_Alltoall exchanges */

typedef struct {

  cs_datatype_t   datatype;          /* associated datatype */
  cs_datatype_t   dest_id_datatype;  /* type of destination id (CS_LNUM_TYPE,
                                        or CS_DATATYPE_NULL) */

  size_t          stride;            /* stride if strided, 0 otherwise */

  size_t          elt_shift;         /* starting byte for element data */
  size_t          comp_size;         /* Composite element size, with padding */

  size_t          send_size;         /* Send buffer element count */
  size_t          recv_size;         /* Receive buffer element count */

  const void     *send_buffer;       /* Send buffer */
  unsigned char  *_send_buffer;      /* Send buffer */

  int            *send_count;        /* Send counts for MPI_Alltoall */
  int            *recv_count;        /* Receive counts for MPI_Alltoall */
  int            *send_displ;        /* Send displs for MPI_Alltoall */
  int            *recv_displ;        /* Receive displs for MPI_Alltoall */

  int            *recv_count_save;   /* Saved (strided) receive counts for
                                        MPI_Alltoall for indexed exchanges */

  MPI_Comm        comm;              /* Associated MPI communicator */
  MPI_Datatype    comp_type;         /* Associated MPI datatype */

  int             n_ranks;           /* Number of ranks associated with
                                        communicator */

} _mpi_all_to_all_caller_t;

/* Base structure for cs_rank_neighbors-based exchanges */

typedef struct {

  cs_rank_neighbors_t   *rn_send;    /* rank neighbors structure for sending */
  cs_rank_neighbors_t   *rn_recv;    /* rank neighbors structure for receiving */

  cs_datatype_t   datatype;          /* associated datatype */
  cs_datatype_t   dest_id_datatype;  /* type of destination id (CS_LNUM_TYPE,
                                        or CS_DATATYPE_NULL) */

  size_t          stride;            /* stride if strided, 0 otherwise */

  size_t          elt_shift;         /* starting byte for element data */
  size_t          comp_size;         /* Composite element size, with padding */

  size_t          send_size;         /* Send buffer element count */
  size_t          recv_size;         /* Receive buffer element count */

  int            *elt_rank_index;    /* Element rank index matching owner
                                        rank_index */

  const void     *send_buffer;       /* Send buffer */
  unsigned char  *_send_buffer;      /* Send buffer */

  int            *send_count;        /* Send counts for MPI_Alltoall */
  int            *recv_count;        /* Receive counts for MPI_Alltoall */
  int            *send_displ;        /* Send displs for MPI_Alltoall */
  int            *recv_displ;        /* Receive displs for MPI_Alltoall */

  int            *recv_count_save;   /* Saved (strided) receive counts for
                                        MPI_Alltoall for indexed exchanges */

  MPI_Comm        comm;              /* Associated MPI communicator */
  MPI_Datatype    comp_type;         /* Associated MPI datatype */

} _hybrid_pex_t;

#endif /* defined(HAVE_MPI) */

/* Structure used to redistribute data */

#if defined(HAVE_MPI)

struct _cs_all_to_all_t {

  cs_lnum_t                  n_elts_src;   /* Number of source elements */
  cs_lnum_t                  n_elts_dest;  /* Number of destination elements
                                              (-1 before metadata available) */
  cs_lnum_t                  n_elts_dest_e; /* Number of destination elements
                                               exchanged, independently of
                                                dest_id (-1 before metadata
                                                available) */

  int                        flags;        /* option flags */

  /* Send metadata */

  const int                 *dest_rank;    /* optional element destination
                                              rank (possibly shared) */
  int                       *_dest_rank;   /* dest_rank if owner, or NULL */

  const cs_lnum_t           *dest_id;      /* optional element destination id
                                              (possibly shared) */
  cs_lnum_t                 *_dest_id;     /* dest_id if owner, or NULL */

  /* Receive metadata */

  cs_lnum_t                 *recv_id;      /* received match for dest_id */

  /* Data needed only for Crystal Router reverse communication */

  cs_lnum_t                 *src_id;       /* received match for dest_id */
  int                       *src_rank;     /* received source rank */

  /* Sub-structures */

  _mpi_all_to_all_caller_t  *dc;       /* Default MPI_Alltoall(v) caller */
  _hybrid_pex_t             *hc;       /* Hybrid PEX */
  cs_crystal_router_t       *cr;       /* associated crystal-router */

  /* MPI data */

  int                        n_ranks;      /* Number of associated ranks */
  MPI_Comm                   comm;         /* associated communicator */

  /* Other metadata */

  cs_all_to_all_type_t       type;         /* Communication protocol */

};

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default all-to-all type */

static cs_all_to_all_type_t _all_to_all_type = CS_ALL_TO_ALL_MPI_DEFAULT;

static cs_rank_neighbors_exchange_t _hybrid_meta_type
  = CS_RANK_NEIGHBORS_CRYSTAL_ROUTER;

#if defined(HAVE_MPI)

/* Call counter and timer: 0: total, 1: metadata comm, 2: data comm */

static size_t              _all_to_all_calls[3] = {0, 0, 0};
static cs_timer_counter_t  _all_to_all_timers[3];

/* Instrumentation */

static int        _n_trace = 0;
static int        _n_trace_max = 0;
static uint64_t  *_all_to_all_trace = NULL;
static FILE      *_all_to_all_trace_bt_log = NULL;

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Common portion of different all-to-all distributor contructors.
 *
 * arguments:
 *   n_elts <-- number of elements
 *   flags  <-- sum of ordering and metadata flag constants
 *   comm   <-- associated MPI communicator
 *
 * returns:
 *   pointer to new all-to-all distributor
 *----------------------------------------------------------------------------*/

static cs_all_to_all_t *
_all_to_all_create_base(size_t    n_elts,
                        int       flags,
                        MPI_Comm  comm)
{
  cs_all_to_all_t *d;

  /* Initialize timers if required */

  if (_all_to_all_calls[0] == 0) {
    int n_timers = sizeof(_all_to_all_timers)/sizeof(_all_to_all_timers[0]);
    for (int i = 0; i < n_timers; i++)
      CS_TIMER_COUNTER_INIT(_all_to_all_timers[i]);

    const char *p = getenv("CS_ALL_TO_ALL_TRACE");
    if (p != NULL) {
      if (atoi(p) > 0) {
        _n_trace_max = atoi(p);
        bft_printf(_("\n cs_all_2_all_trace: %d.\n\n"), _n_trace_max);
        BFT_MALLOC(_all_to_all_trace, _n_trace_max*9, uint64_t);
        _all_to_all_trace_bt_log = fopen("all_to_all_trace_bt.txt", "w");
      }
    }
  }

  /* Check flags */

  if (   (flags & CS_ALL_TO_ALL_USE_DEST_ID)
      && (flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK))
    bft_error(__FILE__, __LINE__, 0,
              "%s: flags may not match both\n"
              "CS_ALL_TO_ALL_USE_DEST_ID and\n"
              "CS_ALL_TO_ALL_ORDER_BY_SRC_RANK.",
              __func__);

  /* Allocate structure */

  BFT_MALLOC(d, 1, cs_all_to_all_t);

  /* Create associated sub-structure */

  d->n_elts_src = n_elts;
  d->n_elts_dest = -1; /* undetermined as yet */
  d->n_elts_dest_e = -1;

  d->flags = flags;

  d->dest_rank = NULL;
  d->_dest_rank = NULL;

  d->dest_id = NULL;
  d->_dest_id = NULL;

  d->recv_id = NULL;

  d->src_id = NULL;
  d->src_rank = NULL;

  d->cr = NULL;
  d->hc = NULL;
  d->dc = NULL;

  d->comm = comm;
  MPI_Comm_size(comm, &(d->n_ranks));

  d->type = _all_to_all_type;

  return d;
}

/*----------------------------------------------------------------------------
 * Compute rank displacement based on count.
 *
 * arguments:
 *   n_ranks <-- number of ranks
 *   count   <-- number of elements per rank (size: n_ranks)
 *   displ   --> element displacement in cumulative array (size: n_ranks)
 *
 * returns:
 *   cumulative count for all ranks
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_compute_displ(int        n_ranks,
               const int  count[],
               int        displ[])
{
  int i;
  cs_lnum_t total_count = 0;

  displ[0] = 0;

  for (i = 0; i < n_ranks; i++)
    displ[i+1] = displ[i] + count[i];

  if (n_ranks > 0)
    total_count = displ[n_ranks-1] + count[n_ranks-1];

  return total_count;
}

/*----------------------------------------------------------------------------
 * First stage of creation for an MPI_Alltoall(v) caller for strided data.
 *
 * parameters:
 *   flags     <-- metadata flags
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _mpi_all_to_all_caller_t *
_alltoall_caller_create_meta(int        flags,
                             MPI_Comm   comm)
{
  _mpi_all_to_all_caller_t *dc;

  /* Allocate structure */

  BFT_MALLOC(dc, 1, _mpi_all_to_all_caller_t);

  dc->datatype = CS_DATATYPE_NULL;

  dc->dest_id_datatype = CS_DATATYPE_NULL;

  if (flags & CS_ALL_TO_ALL_USE_DEST_ID)
    dc->dest_id_datatype = CS_LNUM_TYPE;

  dc->stride = 0;

  dc->send_size = 0;
  dc->recv_size = 0;

  dc->comm = comm;

  MPI_Comm_size(comm, &(dc->n_ranks));

  dc->send_buffer = NULL;
  dc->_send_buffer = NULL;

  BFT_MALLOC(dc->send_count, dc->n_ranks, int);
  BFT_MALLOC(dc->recv_count, dc->n_ranks, int);
  BFT_MALLOC(dc->send_displ, dc->n_ranks + 1, int);
  BFT_MALLOC(dc->recv_displ, dc->n_ranks + 1, int);
  dc->recv_count_save = NULL;

  /* Compute data size and alignment */

  if (dc->dest_id_datatype == CS_LNUM_TYPE)
    dc->elt_shift = sizeof(cs_lnum_t);
  else
    dc->elt_shift = 0;

  size_t align_size = sizeof(cs_lnum_t);

  if (dc->elt_shift % align_size)
    dc->elt_shift += align_size - (dc->elt_shift % align_size);

  dc->comp_size = dc->elt_shift;

  dc->comp_type = MPI_BYTE;

  /* Return pointer to structure */

  return dc;
}

/*----------------------------------------------------------------------------
 * Destroy a MPI_Alltoall(v) caller.
 *
 * parameters:
 *   dc <-> pointer to pointer to MPI_Alltoall(v) caller structure
 *---------------------------------------------------------------------------*/

static void
_alltoall_caller_destroy(_mpi_all_to_all_caller_t **dc)
{
  if (dc != NULL) {
    _mpi_all_to_all_caller_t *_dc = *dc;
    if (_dc->comp_type != MPI_BYTE)
      MPI_Type_free(&(_dc->comp_type));
    BFT_FREE(_dc->_send_buffer);
    BFT_FREE(_dc->recv_count_save);
    BFT_FREE(_dc->recv_displ);
    BFT_FREE(_dc->send_displ);
    BFT_FREE(_dc->recv_count);
    BFT_FREE(_dc->send_count);
    BFT_FREE(*dc);
  }
}

/*----------------------------------------------------------------------------
 * Update all MPI_Alltoall(v) caller metadata.
 *
 * parameters:
 *   dc        <-> distributor caller
 *   datatype  <-- associated datatype
 *   stride    <-- associated stride (0 if indexed)
 *---------------------------------------------------------------------------*/

static void
_alltoall_caller_update_meta(_mpi_all_to_all_caller_t  *dc,
                             cs_datatype_t              datatype,
                             int                        stride)
{
  size_t elt_size = cs_datatype_size[datatype]*stride;
  size_t align_size = sizeof(int);

  /* Free previous associated datatype if needed */

  if (dc->comp_type != MPI_BYTE)
    MPI_Type_free(&(dc->comp_type));

  /* Now update metadata */

  dc->datatype = datatype;
  dc->stride = stride;

  /* Recompute data size and alignment */

  if (dc->dest_id_datatype == CS_LNUM_TYPE) {
    dc->elt_shift = sizeof(cs_lnum_t);
    if (sizeof(cs_lnum_t) > align_size)
      align_size = sizeof(cs_lnum_t);
  }
  else
    dc->elt_shift = 0;

  if (dc->elt_shift % align_size)
    dc->elt_shift += align_size - (dc->elt_shift % align_size);

  if (stride > 0) {
    if (cs_datatype_size[datatype] > align_size)
      align_size = cs_datatype_size[datatype];
    dc->comp_size = dc->elt_shift + elt_size;
  }
  else
    dc->comp_size = dc->elt_shift + cs_datatype_size[datatype];

  if (elt_size % align_size)
    dc->comp_size += align_size - (elt_size % align_size);

  /* Update associated MPI datatype */

  MPI_Type_contiguous(dc->comp_size, MPI_BYTE, &(dc->comp_type));
  MPI_Type_commit(&(dc->comp_type));
}

/*----------------------------------------------------------------------------
 * Save partial metadata before indexed MPI_Alltoall(v) call.
 *
 * parameters:
 *   dc <-> associated MPI_Alltoall(v) caller structure
 *---------------------------------------------------------------------------*/

static  void
_alltoall_caller_save_meta_i(_mpi_all_to_all_caller_t  *dc)
{
  if (dc->recv_count_save == NULL) {
    BFT_MALLOC(dc->recv_count_save, dc->n_ranks, int);
    memcpy(dc->recv_count_save, dc->recv_count, sizeof(int)*(dc->n_ranks));
  }
}

/*----------------------------------------------------------------------------
 * Reset metadate to that used for strided exchanges after indexed
 * MPI_Alltoall(v) call.
 *
 * parameters:
 *   dc        <-> associated MPI_Alltoall(v) caller structure
 *   dest_rank <-- destination rank (in direct mode), used here to reorder
 *                 data in reverse mode.
 *---------------------------------------------------------------------------*/

static  void
_alltoall_caller_reset_meta_i(_mpi_all_to_all_caller_t  *dc,
                              const int                  dest_rank[])
{
  /* Re-count values to send */
  for (int i = 0; i < dc->n_ranks; i++)
    dc->send_count[i] = 0;
  for (size_t j = 0; j < dc->send_size; j++)
    dc->send_count[dest_rank[j]] += 1;

  /* Revert to saved values for receive */
  if (dc->recv_count_save != NULL) {
    memcpy(dc->recv_count, dc->recv_count_save, sizeof(int)*(dc->n_ranks));
    BFT_FREE(dc->recv_count_save);
  }

  /* Recompute associated displacements */
  _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);
  _compute_displ(dc->n_ranks, dc->recv_count, dc->recv_displ);
}

/*----------------------------------------------------------------------------
 * Swap source and destination ranks of all-to-all distributor.
 *
 * parameters:
 *   d <->  pointer to associated all-to-all distributor
 *----------------------------------------------------------------------------*/

static void
_alltoall_caller_swap_src_dest(_mpi_all_to_all_caller_t  *dc)
{
  size_t tmp_size[2] = {dc->send_size, dc->recv_size};
  int *tmp_count = dc->recv_count;
  int *tmp_displ = dc->recv_displ;

  dc->send_size = tmp_size[1];
  dc->recv_size = tmp_size[0];

  dc->recv_count = dc->send_count;
  dc->recv_displ = dc->send_displ;

  dc->send_count = tmp_count;
  dc->send_displ = tmp_displ;
}

/*----------------------------------------------------------------------------
 * Prepare a MPI_Alltoall(v) caller for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   reverse          <-- true if reverse mode
 *   data             <-- element values
 *   dest_id          <-- element destination id, or NULL
 *   recv_id          <-- element receive id (for reverse mode), or NULL
 *   dest_rank        <-- destination rank for each element
 *---------------------------------------------------------------------------*/

static void
_alltoall_caller_prepare_s(_mpi_all_to_all_caller_t  *dc,
                           size_t                     n_elts,
                           int                        stride,
                           cs_datatype_t              datatype,
                           bool                       reverse,
                           const void                *data,
                           const cs_lnum_t           *dest_id,
                           const cs_lnum_t           *recv_id,
                           const int                  dest_rank[])
{
  int i;

  size_t elt_size = cs_datatype_size[datatype]*stride;

  unsigned const char *_data = data;

  assert(data != NULL || (n_elts == 0 || stride == 0));

  _alltoall_caller_update_meta(dc, datatype, stride);

  /* Allocate send buffer */

  if (reverse && recv_id == NULL) {
    BFT_FREE(dc->_send_buffer);
    dc->send_buffer = data;
  }
  else {
    BFT_REALLOC(dc->_send_buffer, dc->send_size*dc->comp_size, unsigned char);
    dc->send_buffer = dc->_send_buffer;
  }

  /* Copy metadata */

  if (dc->dest_id_datatype == CS_LNUM_TYPE) {
    unsigned const char *_dest_id = (unsigned const char *)dest_id;
    const size_t id_size = sizeof(cs_lnum_t);
    for (size_t j = 0; j < n_elts; j++) {
      size_t w_displ = dc->send_displ[dest_rank[j]]*dc->comp_size;
      size_t r_displ = j*id_size;
      dc->send_displ[dest_rank[j]] += 1;
      for (size_t k = 0; k < id_size; k++)
        dc->_send_buffer[w_displ + k] = _dest_id[r_displ + k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }

  /* Copy data; in case of reverse send with destination ids, the
     matching indirection must be applied here. */

  if (!reverse) {
    for (size_t j = 0; j < n_elts; j++) {
      size_t w_displ =   dc->send_displ[dest_rank[j]]*dc->comp_size
                       + dc->elt_shift;
      size_t r_displ = j*elt_size;
      dc->send_displ[dest_rank[j]] += 1;
      for (size_t k = 0; k < elt_size; k++)
        dc->_send_buffer[w_displ + k] = _data[r_displ + k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }
  else if (recv_id != NULL) { /* reverse here */
    for (size_t j = 0; j < dc->send_size; j++) {
      size_t w_displ = dc->comp_size*j + dc->elt_shift;
      size_t r_displ = recv_id[j]*elt_size;
      for (size_t k = 0; k < elt_size; k++)
        dc->_send_buffer[w_displ + k] = _data[r_displ + k];
    }
  }
  else { /* If revert and no recv_id */
    assert(dc->send_buffer == data);
  }
}

/*----------------------------------------------------------------------------
 * Prepare a MPI_Alltoall(v) caller for indexed data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   datatype         <-- type of data considered
 *   reverse          <-- true if reverse mode
 *   src_index        <-- source index
 *   dest_index       <-- destination index
 *   data             <-- element values
 *   recv_id          <-- element receive id (for reverse mode), or NULL
 *   dest_rank        <-- destination rank for each element
 *---------------------------------------------------------------------------*/

static void
_alltoall_caller_prepare_i(_mpi_all_to_all_caller_t  *dc,
                           size_t                     n_elts,
                           cs_datatype_t              datatype,
                           bool                       reverse,
                           const cs_lnum_t            src_index[],
                           const cs_lnum_t            dest_index[],
                           const void                *data,
                           const cs_lnum_t           *recv_id,
                           const int                  dest_rank[])
{
  int i;

  size_t elt_size = cs_datatype_size[datatype];

  unsigned const char *_data = data;

  _alltoall_caller_update_meta(dc, datatype, 0);

  /* Allocate send buffer */

  if (reverse && recv_id == NULL) {
    BFT_FREE(dc->_send_buffer);
    dc->send_buffer = data;
  }
  else {
    size_t n_sub_elts = 0;
    if (reverse == false || recv_id == NULL)
      n_sub_elts = src_index[n_elts];
    else {
      for (size_t j = 0; j < dc->send_size; j++) {
        cs_lnum_t k = recv_id[j];
        n_sub_elts += src_index[k+1] - src_index[k];
      }
    }
    BFT_REALLOC(dc->_send_buffer, n_sub_elts*dc->comp_size, unsigned char);
    dc->send_buffer = dc->_send_buffer;
  }

  /* Update send and receive counts and displacements
     (avoiding additional communication) */

  for (i = 0; i < dc->n_ranks; i++) {
    dc->send_count[i] = 0;
    dc->recv_count[i] = 0;
  }

  if (reverse) {
    if (recv_id != NULL) {
      for (i = 0; i < dc->n_ranks; i++) {
        size_t i_s = dc->send_displ[i];
        size_t i_e = dc->send_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t k = recv_id[j];
          cs_lnum_t n_sub_send = src_index[k+1] - src_index[k];
          dc->send_count[i] += n_sub_send;
        }
      }
    }
    else {
      for (i = 0; i < dc->n_ranks; i++) {
        size_t i_s = dc->send_displ[i];
        size_t i_e = dc->send_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
          dc->send_count[i] += n_sub_send;
        }
      }
    }
  }
  else { /* !reverse */
    for (size_t j = 0; j < dc->send_size; j++) {
      cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
      dc->send_count[dest_rank[j]] += n_sub_send;
    }
  }

  if (reverse) {
    for (size_t j = 0; j < dc->recv_size; j++) {
      cs_lnum_t n_sub_recv = dest_index[j+1] - dest_index[j];
      dc->recv_count[dest_rank[j]] += n_sub_recv;
    }
  }
  else {
    if (recv_id != NULL) {
      for (i = 0; i < dc->n_ranks; i++) {
        size_t i_s = dc->recv_displ[i];
        size_t i_e = dc->recv_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t k = recv_id[j];
          cs_lnum_t n_sub_recv = dest_index[k+1] - dest_index[k];
          dc->recv_count[i] += n_sub_recv;
        }
      }
    }
    else {
      for (i = 0; i < dc->n_ranks; i++) {
        size_t i_s = dc->recv_displ[i];
        size_t i_e = dc->recv_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t n_sub_recv = dest_index[j+1] - dest_index[j];
          dc->recv_count[i] += n_sub_recv;
        }
      }
    }
  }

  _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);
  _compute_displ(dc->n_ranks, dc->recv_count, dc->recv_displ);

  /* Copy data; in case of reverse send with destination ids, the
     matching indirection must be applied here. */

  if (!reverse) {
    for (size_t j = 0; j < n_elts; j++) {
      cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
      size_t w_displ = dc->send_displ[dest_rank[j]]*elt_size;
      size_t r_displ = src_index[j]*elt_size;
      size_t sub_size = elt_size*n_sub_send;
      dc->send_displ[dest_rank[j]] += n_sub_send;
      for (size_t l = 0; l < sub_size; l++)
        dc->_send_buffer[w_displ + l] = _data[r_displ + l];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }
  else if (recv_id != NULL) { /* reverse here */
    size_t w_displ = 0;
    for (size_t j = 0; j < dc->send_size; j++) {
      cs_lnum_t k = recv_id[j];
      cs_lnum_t n_sub_send = src_index[k+1] - src_index[k];
      size_t r_displ = src_index[k]*elt_size;
      size_t sub_size = elt_size*n_sub_send;
      for (size_t l = 0; l < sub_size; l++)
        dc->_send_buffer[w_displ + l] = _data[r_displ + l];
      w_displ += sub_size;
    }

    /* Recompute associated displacements */
    _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);
    _compute_displ(dc->n_ranks, dc->recv_count, dc->recv_displ);


  }
  else { /* If revert and no recv_id */
    assert(dc->send_buffer == data);
  }
}

/*----------------------------------------------------------------------------
 * Exchange metadata for an MPI_Alltoall(v) caller.
 *
 * Order of data from a same source rank is preserved.
 *
 * parameters:
 *   dc        <-> associated MPI_Alltoall(v) caller structure
 *   n_elts    <-- number of elements
 *   dest_rank <-- destination rank for each element
 *---------------------------------------------------------------------------*/

static  void
_alltoall_caller_exchange_meta(_mpi_all_to_all_caller_t  *dc,
                               size_t                     n_elts,
                               const int                  dest_rank[])
{
  /* Count values to send and receive */

  for (int i = 0; i < dc->n_ranks; i++)
    dc->send_count[i] = 0;

  for (size_t j = 0; j < n_elts; j++)
    dc->send_count[dest_rank[j]] += 1;

  dc->send_size = _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);

  /* Exchange counts */

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = 0;
    _all_to_all_trace[_n_trace*9+6] = 0;
    _all_to_all_trace[_n_trace*9+7] = 0;
    _all_to_all_trace[_n_trace*9+8] = 0;
    _n_trace += 1;
  }

  MPI_Alltoall(dc->send_count, 1, MPI_INT,
               dc->recv_count, 1, MPI_INT,
               dc->comm);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_METADATA,
                            &t0, &t1);

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 1;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = 0;
    _all_to_all_trace[_n_trace*9+6] = 0;
    _all_to_all_trace[_n_trace*9+7] = 0;
    _all_to_all_trace[_n_trace*9+8] = 0;
    _n_trace += 1;
  }

  _all_to_all_calls[CS_ALL_TO_ALL_TIME_METADATA] += 1;

  dc->recv_size = _compute_displ(dc->n_ranks, dc->recv_count, dc->recv_displ);
}

/*----------------------------------------------------------------------------
 * Exchange strided data with a MPI_Alltoall(v) caller.
 *
 * parameters:
 *   d         <-> pointer to associated all-to-all distributor
 *   dc        <-> associated MPI_Alltoall(v) caller structure
 *   reverse   <-- true if reverse mode
 *   dest_rank <-- destination rank (in direct mode), used here to reorder
 *                 data in reverse mode.
 *   dest_data <-> destination data buffer, or NULL
 *
 * returns:
 *   pointer to dest_data, or newly allocated buffer
 *---------------------------------------------------------------------------*/

static  void *
_alltoall_caller_exchange_s(cs_all_to_all_t           *d,
                            _mpi_all_to_all_caller_t  *dc,
                            bool                       reverse,
                            const int                  dest_rank[],
                            void                      *dest_data)
{
  size_t elt_size = cs_datatype_size[dc->datatype]*dc->stride;
  unsigned char *_dest_data = dest_data, *_recv_data = dest_data;

  /* Final data buffer */

  if (_dest_data == NULL && dc->recv_size*elt_size > 0)
    BFT_MALLOC(_dest_data, dc->recv_size*elt_size, unsigned char);

  /* Data buffer for MPI exchange (may merge data and metadata) */
  if (   dc->dest_id_datatype == CS_LNUM_TYPE || d->recv_id != NULL
      || reverse)
    BFT_MALLOC(_recv_data, dc->recv_size*dc->comp_size, unsigned char);
  else
    _recv_data = _dest_data;

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    int n_r_send = 0, n_r_recv = 0;
    for (int i = 0; i < dc->n_ranks; i++) {
      if (dc->send_count[i] > 0)
        n_r_send += 1;
      if (dc->recv_count[i] > 0)
        n_r_recv += 1;
    }
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 2;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = dc->send_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = dc->recv_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = n_r_send;
    _all_to_all_trace[_n_trace*9+8] = n_r_recv;
    _n_trace += 1;
  }

  MPI_Alltoallv(dc->send_buffer, dc->send_count, dc->send_displ, dc->comp_type,
                _recv_data, dc->recv_count, dc->recv_displ, dc->comp_type,
                dc->comm);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;

  if (_n_trace < _n_trace_max) {
    int n_r_send = 0, n_r_recv = 0;
    for (int i = 0; i < dc->n_ranks; i++) {
      if (dc->send_count[i] > 0)
        n_r_send += 1;
      if (dc->recv_count[i] > 0)
        n_r_recv += 1;
    }
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = dc->send_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = dc->recv_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = n_r_send;
    _all_to_all_trace[_n_trace*9+8] = n_r_recv;
    _n_trace += 1;
  }

  /* dest id datatype only used for first exchange */

  if (dc->dest_id_datatype == CS_LNUM_TYPE) {
    assert(d->recv_id == NULL);
    BFT_MALLOC(d->recv_id, d->dc->recv_size, cs_lnum_t);
    for (size_t i = 0; i < d->dc->recv_size; i++)
      d->recv_id[i] = -1;
    const unsigned char *sp = _recv_data;
    for (size_t i = 0; i < d->dc->recv_size; i++)
      memcpy(d->recv_id + i,
             sp + d->dc->comp_size*i,
             sizeof(cs_lnum_t));
    dc->dest_id_datatype = CS_DATATYPE_NULL;
    cs_lnum_t dest_id_max = -1;
    for (size_t i = 0; i < d->dc->recv_size; i++) {
      if (d->recv_id[i] > dest_id_max)
        dest_id_max = d->recv_id[i];
      d->n_elts_dest = dest_id_max + 1;
    }
    d->n_elts_dest_e = d->dc->recv_size;
  }

  /* Now handle main data buffer (reverse implies reordering data) */

  if (_dest_data != _recv_data) {
    const unsigned char *sp = _recv_data + d->dc->elt_shift;
    if (d->recv_id != NULL && !reverse) {
      for (size_t i = 0; i < d->dc->recv_size; i++) {
        size_t w_displ = d->recv_id[i]*elt_size;
        size_t r_displ = d->dc->comp_size*i;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
      }
    }
    else if (reverse) {
      for (int i = 0; i < dc->n_ranks; i++)
        dc->recv_count[i] = 0;
      for (size_t i = 0; i < d->dc->recv_size; i++) {
        int rank_id = dest_rank[i];
        size_t w_displ = i*elt_size;
        size_t r_displ = (  dc->recv_displ[rank_id]
                          + dc->recv_count[rank_id])*dc->comp_size;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
        dc->recv_count[rank_id] += 1;
      }
    }
    else {
      for (size_t i = 0; i < d->dc->recv_size; i++) {
        size_t w_displ = i*elt_size;
        size_t r_displ = dc->comp_size*i;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
      }
    }
    BFT_FREE(_recv_data);
  }

  return _dest_data;
}

/*----------------------------------------------------------------------------
 * Exchange indexed data with a MPI_Alltoall(v) caller.
 *
 * parameters:
 *   d          <-> pointer to associated all-to-all distributor
 *   dc         <-> associated MPI_Alltoall(v) caller structure
 *   reverse    <-- true if reverse mode
 *   dest_index <-- destination index
 *   dest_rank  <-- destination rank (in direct mode), used here to reorder
 *                  data in reverse mode.
 *   dest_data  <-> destination data buffer, or NULL
 *
 * returns:
 *   pointer to dest_data, or newly allocated buffer
 *---------------------------------------------------------------------------*/

static  void *
_alltoall_caller_exchange_i(cs_all_to_all_t           *d,
                            _mpi_all_to_all_caller_t  *dc,
                            bool                       reverse,
                            const cs_lnum_t            dest_index[],
                            const int                  dest_rank[],
                            void                      *dest_data)
{
  size_t elt_size = cs_datatype_size[dc->datatype];
  unsigned char *_dest_data = dest_data, *_recv_data = dest_data;

  /* Final data buffer */

  size_t n_elts_dest = (reverse) ? d->n_elts_src : d->n_elts_dest;

  size_t _n_dest_sub = dest_index[n_elts_dest];
  if (_dest_data == NULL && _n_dest_sub*elt_size > 0)
    BFT_MALLOC(_dest_data, _n_dest_sub*elt_size, unsigned char);

  /* Data buffer for MPI exchange (may merge data and metadata) */

  if (d->recv_id != NULL || reverse) {
    size_t _n_dest_buf = dc->recv_displ[dc->n_ranks] * dc->comp_size;
    BFT_MALLOC(_recv_data, _n_dest_buf, unsigned char);
  }
  else
    _recv_data = _dest_data;

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    int n_r_send = 0, n_r_recv = 0;
    for (int i = 0; i < dc->n_ranks; i++) {
      if (dc->send_count[i] > 0)
        n_r_send += 1;
      if (dc->recv_count[i] > 0)
        n_r_recv += 1;
    }
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 2;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = dc->send_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = dc->recv_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = n_r_send;
    _all_to_all_trace[_n_trace*9+8] = n_r_recv;
    _n_trace += 1;
  }

  MPI_Alltoallv(dc->send_buffer, dc->send_count, dc->send_displ, dc->comp_type,
                _recv_data, dc->recv_count, dc->recv_displ, dc->comp_type,
                dc->comm);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;

  if (_n_trace < _n_trace_max) {
    int n_r_send = 0, n_r_recv = 0;
    for (int i = 0; i < dc->n_ranks; i++) {
      if (dc->send_count[i] > 0)
        n_r_send += 1;
      if (dc->recv_count[i] > 0)
        n_r_recv += 1;
    }
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = dc->send_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = dc->recv_size*dc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = n_r_send;
    _all_to_all_trace[_n_trace*9+8] = n_r_recv;
    _n_trace += 1;
  }

  /* Handle main data buffer (reverse implies reordering data) */

  if (_dest_data != _recv_data) {
    const unsigned char *sp = _recv_data;
    if (d->recv_id != NULL && !reverse) {
      size_t r_displ = 0;
      for (size_t i = 0; i < dc->recv_size; i++) {
        cs_lnum_t k = d->recv_id[i];
        cs_lnum_t n_sub_recv = dest_index[k+1] - dest_index[k];
        size_t w_displ = dest_index[k]*elt_size;
        size_t sub_size = n_sub_recv*elt_size;
        for (size_t l = 0; l < sub_size; l++)
          _dest_data[w_displ + l] = sp[r_displ + l];
        r_displ += sub_size;
      }
    }
    else if (reverse) {
      for (int i = 0; i < dc->n_ranks; i++)
        dc->recv_count[i] = 0;
      for (size_t i = 0; i < d->dc->recv_size; i++) {
        int rank_id = dest_rank[i];
        size_t w_displ = dest_index[i]*elt_size;
        size_t r_displ = (  dc->recv_displ[rank_id]
                          + dc->recv_count[rank_id])*elt_size;
        cs_lnum_t n_sub_recv = dest_index[i+1] - dest_index[i];
        size_t sub_size = n_sub_recv*elt_size;
        for (size_t l = 0; l < sub_size; l++)
          _dest_data[w_displ + l] = sp[r_displ + l];
        dc->recv_count[rank_id] += n_sub_recv;
      }
    }
    else {
      assert(_dest_data == _recv_data);
    }
    BFT_FREE(_recv_data);
  }

  return _dest_data;
}

/*----------------------------------------------------------------------------
 * First stage of creation for a hybrid caller for strided data.
 *
 * parameters:
 *   flags     <-- metadata flags
 *   comm      <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _hybrid_pex_t *
_hybrid_pex_create_meta(int        flags,
                        MPI_Comm   comm)
{
  _hybrid_pex_t *hc;

  /* Allocate structure */

  BFT_MALLOC(hc, 1, _hybrid_pex_t);

  hc->rn_send = NULL;
  hc->rn_recv = NULL;

  hc->datatype = CS_DATATYPE_NULL;

  hc->dest_id_datatype = CS_DATATYPE_NULL;

  if (flags & CS_ALL_TO_ALL_USE_DEST_ID)
    hc->dest_id_datatype = CS_LNUM_TYPE;

  hc->stride = 0;

  hc->send_size = 0;
  hc->recv_size = 0;

  hc->comm = comm;

  hc->elt_rank_index = NULL;

  hc->send_buffer = NULL;
  hc->_send_buffer = NULL;

  hc->send_count = NULL;
  hc->recv_count = NULL;
  hc->send_displ = NULL;
  hc->recv_displ = NULL;
  hc->recv_count_save = NULL;

  /* Compute data size and alignment */

  if (hc->dest_id_datatype == CS_LNUM_TYPE)
    hc->elt_shift = sizeof(cs_lnum_t);
  else
    hc->elt_shift = 0;

  size_t align_size = sizeof(cs_lnum_t);

  if (hc->elt_shift % align_size)
    hc->elt_shift += align_size - (hc->elt_shift % align_size);

  hc->comp_size = hc->elt_shift;

  hc->comp_type = MPI_BYTE;

  /* Return pointer to structure */

  return hc;
}

/*----------------------------------------------------------------------------
 * Destroy a hybrid caller.
 *
 * parameters:
 *   hc <-> pointer to pointer to MPI_Alltoall(v) caller structure
 *---------------------------------------------------------------------------*/

static void
_hybrid_pex_destroy(_hybrid_pex_t **hc)
{
  if (hc != NULL) {
    _hybrid_pex_t *_hc = *hc;
    if (_hc->comp_type != MPI_BYTE)
      MPI_Type_free(&(_hc->comp_type));
    BFT_FREE(_hc->elt_rank_index);
    BFT_FREE(_hc->_send_buffer);
    BFT_FREE(_hc->recv_count_save);
    BFT_FREE(_hc->recv_displ);
    BFT_FREE(_hc->send_displ);
    BFT_FREE(_hc->recv_count);
    BFT_FREE(_hc->send_count);

    cs_rank_neighbors_destroy(&(_hc->rn_send));
    cs_rank_neighbors_destroy(&(_hc->rn_recv));

    BFT_FREE(*hc);
  }
}

/*----------------------------------------------------------------------------
 * Update all hybrid caller metadata.
 *
 * parameters:
 *   hc        <-> distributor caller
 *   datatype  <-- associated datatype
 *   stride    <-- associated stride (0 if indexed)
 *---------------------------------------------------------------------------*/

static void
_hybrid_pex_update_meta(_hybrid_pex_t  *hc,
                        cs_datatype_t   datatype,
                        int             stride)
{
  size_t elt_size = cs_datatype_size[datatype]*stride;
  size_t align_size = sizeof(int);

  /* Free previous associated datatype if needed */

  if (hc->comp_type != MPI_BYTE)
    MPI_Type_free(&(hc->comp_type));

  /* Now update metadata */

  hc->datatype = datatype;
  hc->stride = stride;

  /* Recompute data size and alignment */

  if (hc->dest_id_datatype == CS_LNUM_TYPE) {
    hc->elt_shift = sizeof(cs_lnum_t);
    if (sizeof(cs_lnum_t) > align_size)
      align_size = sizeof(cs_lnum_t);
  }
  else
    hc->elt_shift = 0;

  if (hc->elt_shift % align_size)
    hc->elt_shift += align_size - (hc->elt_shift % align_size);

  if (stride > 0) {
    if (cs_datatype_size[datatype] > align_size)
      align_size = cs_datatype_size[datatype];
    hc->comp_size = hc->elt_shift + elt_size;
  }
  else
    hc->comp_size = hc->elt_shift + cs_datatype_size[datatype];

  if (elt_size % align_size)
    hc->comp_size += align_size - (elt_size % align_size);

  /* Update associated MPI datatype */

  MPI_Type_contiguous(hc->comp_size, MPI_BYTE, &(hc->comp_type));
  MPI_Type_commit(&(hc->comp_type));
}

/*----------------------------------------------------------------------------
 * Save partial metadata before indexed hybrid call.
 *
 * parameters:
 *   hc  <-> associated hybrid caller structure
 *---------------------------------------------------------------------------*/

static  void
_hybrid_pex_save_meta_i(_hybrid_pex_t  *hc)
{
  if (hc->recv_count_save == NULL) {
    int n_n_ranks = hc->rn_recv->size;
    BFT_MALLOC(hc->recv_count_save, n_n_ranks, int);
    memcpy(hc->recv_count_save, hc->recv_count, sizeof(int)*(n_n_ranks));
  }
}

/*----------------------------------------------------------------------------
 * Reset metadate to that used for strided exchanges after indexed
 * hybrid call.
 *
 * parameters:
 *   hc              <-> associated hybrid caller structure
 *---------------------------------------------------------------------------*/

static  void
_hybrid_pex_reset_meta_i(_hybrid_pex_t  *hc)
{
  int n_s_ranks = hc->rn_send->size;
  int n_r_ranks = hc->rn_recv->size;

  const int  *dest_rank_index = hc->elt_rank_index;

  /* Re-count values to send */
  for (int i = 0; i < n_s_ranks; i++)
    hc->send_count[i] = 0;
  for (size_t j = 0; j < hc->send_size; j++)
    hc->send_count[dest_rank_index[j]] += 1;

  /* Revert to saved values for receive */
  if (hc->recv_count_save != NULL) {
    memcpy(hc->recv_count, hc->recv_count_save, sizeof(int)*(n_r_ranks));
    BFT_FREE(hc->recv_count_save);
  }

  /* Recompute associated displacements */
  _compute_displ(n_s_ranks, hc->send_count, hc->send_displ);
  _compute_displ(n_r_ranks, hc->recv_count, hc->recv_displ);
}

/*----------------------------------------------------------------------------
 * Swap source and destination ranks of all-to-all distributor.
 *
 * parameters:
 *   d <->  pointer to associated all-to-all distributor
 *----------------------------------------------------------------------------*/

static void
_hybrid_pex_swap_src_dest(_hybrid_pex_t  *hc)
{
  size_t tmp_size[2] = {hc->send_size, hc->recv_size};
  int *tmp_count = hc->recv_count;
  int *tmp_displ = hc->recv_displ;
  cs_rank_neighbors_t *tmp_rn = hc->rn_send;

  hc->rn_send = hc->rn_recv;
  hc->rn_recv = tmp_rn;

  hc->send_size = tmp_size[1];
  hc->recv_size = tmp_size[0];

  hc->recv_count = hc->send_count;
  hc->recv_displ = hc->send_displ;

  hc->send_count = tmp_count;
  hc->send_displ = tmp_displ;
}

/*----------------------------------------------------------------------------
 * Prepare a hybrid caller for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   reverse          <-- true if reverse mode
 *   data             <-- element values
 *   dest_id          <-- element destination id, or NULL
 *   recv_id          <-- element receive id (for reverse mode), or NULL
 *---------------------------------------------------------------------------*/

static void
_hybrid_pex_prepare_s(_hybrid_pex_t    *hc,
                      size_t            n_elts,
                      int               stride,
                      cs_datatype_t     datatype,
                      bool              reverse,
                      const void       *data,
                      const cs_lnum_t  *dest_id,
                      const cs_lnum_t  *recv_id)
{
  int i;

  size_t elt_size = cs_datatype_size[datatype]*stride;

  unsigned const char *_data = data;

  assert(data != NULL || (n_elts == 0 || stride == 0));

  _hybrid_pex_update_meta(hc, datatype, stride);

  /* Allocate send buffer */

  if (reverse && recv_id == NULL) {
    BFT_FREE(hc->_send_buffer);
    hc->send_buffer = data;
  }
  else {
    BFT_REALLOC(hc->_send_buffer, hc->send_size*hc->comp_size, unsigned char);
    hc->send_buffer = hc->_send_buffer;
  }

  const int  *dest_rank_index = hc->elt_rank_index;

  /* Copy metadata */

  if (hc->dest_id_datatype == CS_LNUM_TYPE) {
    unsigned const char *_dest_id = (unsigned const char *)dest_id;
    const size_t id_size = sizeof(cs_lnum_t);
    for (size_t j = 0; j < n_elts; j++) {
      size_t w_displ = hc->send_displ[dest_rank_index[j]]*hc->comp_size;
      size_t r_displ = j*id_size;
      hc->send_displ[dest_rank_index[j]] += 1;
      for (size_t k = 0; k < id_size; k++)
        hc->_send_buffer[w_displ + k] = _dest_id[r_displ + k];
    }
    /* Reset send_displ */
    int n_s_ranks = hc->rn_send->size;
    for (i = 0; i < n_s_ranks; i++)
      hc->send_displ[i] -= hc->send_count[i];
  }

  /* Copy data; in case of reverse send with destination ids, the
     matching indirection must be applied here. */

  if (!reverse) {
    for (size_t j = 0; j < n_elts; j++) {
      size_t w_displ =   hc->send_displ[dest_rank_index[j]]*hc->comp_size
                       + hc->elt_shift;
      size_t r_displ = j*elt_size;
      hc->send_displ[dest_rank_index[j]] += 1;
      for (size_t k = 0; k < elt_size; k++)
        hc->_send_buffer[w_displ + k] = _data[r_displ + k];
    }
    /* Reset send_displ */
    int n_s_ranks = hc->rn_send->size;
    for (i = 0; i < n_s_ranks; i++)
      hc->send_displ[i] -= hc->send_count[i];
  }
  else if (recv_id != NULL) { /* reverse here */
    for (size_t j = 0; j < hc->send_size; j++) {
      size_t w_displ = hc->comp_size*j + hc->elt_shift;
      size_t r_displ = recv_id[j]*elt_size;
      for (size_t k = 0; k < elt_size; k++)
        hc->_send_buffer[w_displ + k] = _data[r_displ + k];
    }
  }
  else { /* If revert and no recv_id */
    assert(hc->send_buffer == data);
  }
}

/*----------------------------------------------------------------------------
 * Prepare a hybrid caller for indexed data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   datatype         <-- type of data considered
 *   reverse          <-- true if reverse mode
 *   src_index        <-- source index
 *   dest_index       <-- destination index
 *   data             <-- element values
 *   recv_id          <-- element receive id (for reverse mode), or NULL
 *---------------------------------------------------------------------------*/

static void
_hybrid_pex_prepare_i(_hybrid_pex_t      *hc,
                      size_t              n_elts,
                      cs_datatype_t       datatype,
                      bool                reverse,
                      const cs_lnum_t     src_index[],
                      const cs_lnum_t     dest_index[],
                      const void         *data,
                      const cs_lnum_t    *recv_id)
{
  int i;

  size_t elt_size = cs_datatype_size[datatype];

  unsigned const char *_data = data;

  _hybrid_pex_update_meta(hc, datatype, 0);

  /* Allocate send buffer */

  if (reverse && recv_id == NULL) {
    BFT_FREE(hc->_send_buffer);
    hc->send_buffer = data;
  }
  else {
    size_t n_sub_elts = 0;
    if (reverse == false || recv_id == NULL)
      n_sub_elts = src_index[n_elts];
    else {
      for (size_t j = 0; j < hc->send_size; j++) {
        cs_lnum_t k = recv_id[j];
        n_sub_elts += src_index[k+1] - src_index[k];
      }
    }
    BFT_REALLOC(hc->_send_buffer, n_sub_elts*hc->comp_size, unsigned char);
    hc->send_buffer = hc->_send_buffer;
  }

  /* Update send and receive counts and displacements
     (avoiding additional communication) */

  const int n_s_ranks = hc->rn_send->size;
  const int n_r_ranks = hc->rn_recv->size;

  for (i = 0; i < n_s_ranks; i++)
    hc->send_count[i] = 0;

  for (i = 0; i < n_r_ranks; i++)
    hc->recv_count[i] = 0;

  const int  *dest_rank_index = hc->elt_rank_index;

  if (reverse) {
    if (recv_id != NULL) {
      for (i = 0; i < n_s_ranks; i++) {
        size_t i_s = hc->send_displ[i];
        size_t i_e = hc->send_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t k = recv_id[j];
          cs_lnum_t n_sub_send = src_index[k+1] - src_index[k];
          hc->send_count[i] += n_sub_send;
        }
      }
    }
    else {
      for (i = 0; i < n_s_ranks; i++) {
        size_t i_s = hc->send_displ[i];
        size_t i_e = hc->send_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
          hc->send_count[i] += n_sub_send;
        }
      }
    }
  }
  else { /* !reverse */
    for (size_t j = 0; j < hc->send_size; j++) {
      cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
      hc->send_count[dest_rank_index[j]] += n_sub_send;
    }
  }

  if (reverse) {
    for (size_t j = 0; j < hc->recv_size; j++) {
      cs_lnum_t n_sub_recv = dest_index[j+1] - dest_index[j];
      hc->recv_count[dest_rank_index[j]] += n_sub_recv;
    }
  }
  else {
    if (recv_id != NULL) {
      for (i = 0; i < n_r_ranks; i++) {
        size_t i_s = hc->recv_displ[i];
        size_t i_e = hc->recv_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t k = recv_id[j];
          cs_lnum_t n_sub_recv = dest_index[k+1] - dest_index[k];
          hc->recv_count[i] += n_sub_recv;
        }
      }
    }
    else {
      for (i = 0; i < n_r_ranks; i++) {
        size_t i_s = hc->recv_displ[i];
        size_t i_e = hc->recv_displ[i+1];
        for (size_t j = i_s; j < i_e; j++) {
          cs_lnum_t n_sub_recv = dest_index[j+1] - dest_index[j];
          hc->recv_count[i] += n_sub_recv;
        }
      }
    }
  }

  _compute_displ(n_s_ranks, hc->send_count, hc->send_displ);
  _compute_displ(n_r_ranks, hc->recv_count, hc->recv_displ);

  /* Copy data; in case of reverse send with destination ids, the
     matching indirection must be applied here. */

  if (!reverse) {
    for (size_t j = 0; j < n_elts; j++) {
      cs_lnum_t n_sub_send = src_index[j+1] - src_index[j];
      size_t w_displ = hc->send_displ[dest_rank_index[j]]*elt_size;
      size_t r_displ = src_index[j]*elt_size;
      size_t sub_size = elt_size*n_sub_send;
      hc->send_displ[dest_rank_index[j]] += n_sub_send;
      for (size_t l = 0; l < sub_size; l++)
        hc->_send_buffer[w_displ + l] = _data[r_displ + l];
    }
    /* Reset send_displ */
    for (i = 0; i < n_s_ranks; i++)
      hc->send_displ[i] -= hc->send_count[i];
  }
  else if (recv_id != NULL) { /* reverse here */
    size_t w_displ = 0;
    for (size_t j = 0; j < hc->send_size; j++) {
      cs_lnum_t k = recv_id[j];
      cs_lnum_t n_sub_send = src_index[k+1] - src_index[k];
      size_t r_displ = src_index[k]*elt_size;
      size_t sub_size = elt_size*n_sub_send;
      for (size_t l = 0; l < sub_size; l++)
        hc->_send_buffer[w_displ + l] = _data[r_displ + l];
      w_displ += sub_size;
    }

    /* Recompute associated displacements */
    _compute_displ(n_s_ranks, hc->send_count, hc->send_displ);
    _compute_displ(n_r_ranks, hc->recv_count, hc->recv_displ);


  }
  else { /* If revert and no recv_id */
    assert(hc->send_buffer == data);
  }
}

/*----------------------------------------------------------------------------
 * Exchange metadata for an hybrid caller.
 *
 * Order of data from a same source rank is preserved.
 *
 * parameters:
 *   hc        <-> associated hybrid caller structure
 *   n_elts    <-- number of elements
 *   dest_rank <-- destination rank for each element
 *---------------------------------------------------------------------------*/

static  void
_hybrid_pex_exchange_meta(_hybrid_pex_t  *hc,
                          size_t          n_elts,
                          const int       dest_rank[])
{
  if (hc->rn_send == NULL) {
    hc->rn_send = cs_rank_neighbors_create(n_elts, dest_rank);

    BFT_MALLOC(hc->elt_rank_index, n_elts, int);
    BFT_MALLOC(hc->send_count, hc->rn_send->size, int);
    BFT_MALLOC(hc->send_displ, hc->rn_send->size + 1, int);

    cs_rank_neighbors_to_index(hc->rn_send,
                               n_elts,
                               dest_rank,
                               hc->elt_rank_index);

  }

  /* Count values to send and receive */

  cs_rank_neighbors_count(hc->rn_send,
                          n_elts,
                          hc->elt_rank_index,
                          hc->send_count);

  hc->send_size = _compute_displ(hc->rn_send->size,
                                 hc->send_count,
                                 hc->send_displ);

  /* Exchange counts */

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = 0;
    _all_to_all_trace[_n_trace*9+6] = 0;
    _all_to_all_trace[_n_trace*9+7] = 0;
    _all_to_all_trace[_n_trace*9+8] = 0;
    _n_trace += 1;
  }

  assert(hc->rn_recv == NULL);

  if (hc->rn_recv != NULL) {
    cs_rank_neighbors_destroy(&(hc->rn_recv));
    BFT_FREE(hc->recv_count);
    BFT_FREE(hc->recv_displ);
  }

  cs_rank_neighbors_sync_count_m(hc->rn_send,
                                 &(hc->rn_recv),
                                 hc->send_count,
                                 &(hc->recv_count),
                                 _hybrid_meta_type,
                                 hc->comm);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_METADATA,
                            &t0, &t1);

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 1;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = 0;
    _all_to_all_trace[_n_trace*9+6] = 0;
    _all_to_all_trace[_n_trace*9+7] = 0;
    _all_to_all_trace[_n_trace*9+8] = 0;
    _n_trace += 1;
  }

  _all_to_all_calls[CS_ALL_TO_ALL_TIME_METADATA] += 1;

  BFT_MALLOC(hc->recv_displ, hc->rn_recv->size + 1, int);

  hc->recv_size = _compute_displ(hc->rn_recv->size,
                                 hc->recv_count,
                                 hc->recv_displ);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange indexed data with a hybrid caller.
 *
 * \param[in]    hc          associated hybrid caller structure
 * \param[in]    sendbuf     starting address of send buffer.
 * \param[out]   recvbuf     starting address of receive buffer.
 */
/*----------------------------------------------------------------------------*/

static  void
_hybrid_alltoallv(_hybrid_pex_t   *hc,
                  const void      *sendbuf,
                  void            *recvbuf)
{
  /* Currently available: MPI_Alltoallv */

  if (true) {
    int n_ranks;
    MPI_Comm_size(hc->comm, &n_ranks);

    int  *send_count, *send_displ, *recv_count, *recv_displ;

    BFT_MALLOC(send_count, n_ranks, int);
    BFT_MALLOC(send_displ, n_ranks+1, int);
    BFT_MALLOC(recv_count, n_ranks, int);
    BFT_MALLOC(recv_displ, n_ranks+1, int);

    const int n_s_ranks = hc->rn_send->size;
    const int n_r_ranks = hc->rn_recv->size;
    const int *s_rank = hc->rn_send->rank;
    const int *r_rank = hc->rn_recv->rank;

    for (int i = 0; i < n_ranks; i++) {
      send_count[i] = 0;
      recv_count[i] = 0;
    }
    for (int i = 0; i < n_s_ranks; i++)
      send_count[s_rank[i]] = hc->send_displ[i+1] - hc->send_displ[i];
    for (int i = 0; i < n_r_ranks; i++)
      recv_count[r_rank[i]] = hc->recv_displ[i+1] - hc->recv_displ[i];

    _compute_displ(n_ranks, send_count, send_displ);
    _compute_displ(n_ranks, recv_count, recv_displ);

    MPI_Alltoallv(sendbuf, send_count, send_displ, hc->comp_type,
                  recvbuf, recv_count, recv_displ, hc->comp_type,
                  hc->comm);

    BFT_FREE(recv_displ);
    BFT_FREE(recv_count);
    BFT_FREE(send_displ);
    BFT_FREE(send_count);
  }
}

/*----------------------------------------------------------------------------
 * Exchange strided data with a hybrid caller.
 *
 * parameters:
 *   d               <-> pointer to associated all-to-all distributor
 *   hc              <-> associated hybrid caller structure
 *   reverse         <-- true if reverse mode
 *   dest_data       <-> destination data buffer, or NULL
 *
 * returns:
 *   pointer to dest_data, or newly allocated buffer
 *---------------------------------------------------------------------------*/

static  void *
_hybrid_pex_exchange_s(cs_all_to_all_t  *d,
                       _hybrid_pex_t    *hc,
                       bool              reverse,
                       void             *dest_data)
{
  size_t elt_size = cs_datatype_size[hc->datatype]*hc->stride;
  unsigned char *_dest_data = dest_data, *_recv_data = dest_data;

  /* Final data buffer */

  if (_dest_data == NULL && hc->recv_size*elt_size > 0)
    BFT_MALLOC(_dest_data, hc->recv_size*elt_size, unsigned char);

  /* Data buffer for MPI exchange (may merge data and metadata) */
  if (   hc->dest_id_datatype == CS_LNUM_TYPE || d->recv_id != NULL
      || reverse)
    BFT_MALLOC(_recv_data, hc->recv_size*hc->comp_size, unsigned char);
  else
    _recv_data = _dest_data;

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 2;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = hc->send_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = hc->recv_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = hc->rn_send->size;
    _all_to_all_trace[_n_trace*9+8] = hc->rn_recv->size;
    _n_trace += 1;
  }

  _hybrid_alltoallv(hc, hc->send_buffer, _recv_data);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = hc->send_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = hc->recv_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = hc->rn_send->size;
    _all_to_all_trace[_n_trace*9+8] = hc->rn_recv->size;
    _n_trace += 1;
  }

  /* dest id datatype only used for first exchange */

  if (hc->dest_id_datatype == CS_LNUM_TYPE) {
    assert(d->recv_id == NULL);
    BFT_MALLOC(d->recv_id, d->hc->recv_size, cs_lnum_t);
    for (size_t i = 0; i < d->hc->recv_size; i++)
      d->recv_id[i] = -1;
    const unsigned char *sp = _recv_data;
    for (size_t i = 0; i < d->hc->recv_size; i++)
      memcpy(d->recv_id + i,
             sp + d->hc->comp_size*i,
             sizeof(cs_lnum_t));
    hc->dest_id_datatype = CS_DATATYPE_NULL;
    cs_lnum_t dest_id_max = -1;
    for (size_t i = 0; i < d->hc->recv_size; i++) {
      if (d->recv_id[i] > dest_id_max)
        dest_id_max = d->recv_id[i];
      d->n_elts_dest = dest_id_max + 1;
    }
    d->n_elts_dest_e = d->hc->recv_size;
  }

  /* Now handle main data buffer (reverse implies reordering data) */

  if (_dest_data != _recv_data) {
    const unsigned char *sp = _recv_data + d->hc->elt_shift;
    if (d->recv_id != NULL && !reverse) {
      for (size_t i = 0; i < d->hc->recv_size; i++) {
        size_t w_displ = d->recv_id[i]*elt_size;
        size_t r_displ = d->hc->comp_size*i;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
      }
    }
    else if (reverse) {
      const int n_r_ranks = hc->rn_recv->size;
      const int *dest_rank_index = hc->elt_rank_index;
      for (int i = 0; i < n_r_ranks; i++)
        hc->recv_count[i] = 0;
      for (size_t i = 0; i < d->hc->recv_size; i++) {
        int rank_idx = dest_rank_index[i];
        size_t w_displ = i*elt_size;
        size_t r_displ = (  hc->recv_displ[rank_idx]
                          + hc->recv_count[rank_idx])*hc->comp_size;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
        hc->recv_count[rank_idx] += 1;
      }
    }
    else {
      for (size_t i = 0; i < d->hc->recv_size; i++) {
        size_t w_displ = i*elt_size;
        size_t r_displ = hc->comp_size*i;
        for (size_t j = 0; j < elt_size; j++)
          _dest_data[w_displ + j] = sp[r_displ + j];
      }
    }
    BFT_FREE(_recv_data);
  }

  return _dest_data;
}

/*----------------------------------------------------------------------------
 * Exchange indexed data with a hybrid caller.
 *
 * parameters:
 *   d                <-> pointer to associated all-to-all distributor
 *   hc               <-> associated hybrid caller structure
 *   reverse          <-- true if reverse mode
 *   dest_index       <-- destination index
 *   dest_data        <-> destination data buffer, or NULL
 *
 * returns:
 *   pointer to dest_data, or newly allocated buffer
 *---------------------------------------------------------------------------*/

static  void *
_hybrid_pex_exchange_i(cs_all_to_all_t    *d,
                       _hybrid_pex_t      *hc,
                       bool                reverse,
                       const cs_lnum_t     dest_index[],
                       void               *dest_data)
{
  size_t elt_size = cs_datatype_size[hc->datatype];
  unsigned char *_dest_data = dest_data, *_recv_data = dest_data;

  /* Final data buffer */

  size_t n_elts_dest = (reverse) ? d->n_elts_src : d->n_elts_dest;

  size_t _n_dest_sub = dest_index[n_elts_dest];
  if (_dest_data == NULL && _n_dest_sub*elt_size > 0)
    BFT_MALLOC(_dest_data, _n_dest_sub*elt_size, unsigned char);

  /* Data buffer for MPI exchange (may merge data and metadata) */

  if (d->recv_id != NULL || reverse) {
    int n_r_ranks = hc->rn_recv->size;
    size_t _n_dest_buf = hc->recv_displ[n_r_ranks] * hc->comp_size;
    BFT_MALLOC(_recv_data, _n_dest_buf, unsigned char);
  }
  else
    _recv_data = _dest_data;

  cs_timer_t t0 = cs_timer_time();

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t0.sec*1e5 + t0.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] = 0;
    _all_to_all_trace[_n_trace*9+2] = 2;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = hc->send_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = hc->recv_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = hc->rn_send->size;
    _all_to_all_trace[_n_trace*9+8] = hc->rn_recv->size;
    _n_trace += 1;
  }

  _hybrid_alltoallv(hc, hc->send_buffer, _recv_data);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;

  if (_n_trace < _n_trace_max) {
    /* Time to 1-5 s */
    _all_to_all_trace[_n_trace*9] = t1.sec*1e5 + t1.nsec/1e4;
    _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                      - _all_to_all_trace[(_n_trace-1)*9];
    _all_to_all_trace[_n_trace*9+2] = 0;
    _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
    _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
    _all_to_all_trace[_n_trace*9+5] = hc->send_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+6] = hc->recv_size*hc->comp_size;
    _all_to_all_trace[_n_trace*9+7] = hc->rn_send->size;
    _all_to_all_trace[_n_trace*9+8] = hc->rn_recv->size;
    _n_trace += 1;
  }

  /* Handle main data buffer (reverse implies reordering data) */

  if (_dest_data != _recv_data) {
    const unsigned char *sp = _recv_data;
    if (d->recv_id != NULL && !reverse) {
      size_t r_displ = 0;
      for (size_t i = 0; i < hc->recv_size; i++) {
        cs_lnum_t k = d->recv_id[i];
        cs_lnum_t n_sub_recv = dest_index[k+1] - dest_index[k];
        size_t w_displ = dest_index[k]*elt_size;
        size_t sub_size = n_sub_recv*elt_size;
        for (size_t l = 0; l < sub_size; l++)
          _dest_data[w_displ + l] = sp[r_displ + l];
        r_displ += sub_size;
      }
    }
    else if (reverse) {
      const int n_r_ranks = hc->rn_recv->size;
      const int *dest_rank_index = hc->elt_rank_index;
      for (int i = 0; i < n_r_ranks; i++)
        hc->recv_count[i] = 0;
      for (size_t i = 0; i < d->hc->recv_size; i++) {
        int rank_idx = dest_rank_index[i];
        size_t w_displ = dest_index[i]*elt_size;
        size_t r_displ = (  hc->recv_displ[rank_idx]
                          + hc->recv_count[rank_idx])*elt_size;
        cs_lnum_t n_sub_recv = dest_index[i+1] - dest_index[i];
        size_t sub_size = n_sub_recv*elt_size;
        for (size_t l = 0; l < sub_size; l++)
          _dest_data[w_displ + l] = sp[r_displ + l];
        hc->recv_count[rank_idx] += n_sub_recv;
      }
    }
    else {
      assert(_dest_data == _recv_data);
    }
    BFT_FREE(_recv_data);
  }

  return _dest_data;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute receive id by ordering based on source rank.
 *
 * The caller is responsible for freeing the returned array.
 *
 * As communication schemes (such as Crystal Router) should lead to
 * elements being grouped by rank, this is used here to optimize sorting.
 * Sorting should also be stable as a result.
 *
 * \param[in]  d         pointer to all-to-all distributor
 * \param[in]  src_rank  associated source rank (size: d->n_elts_dest)
 */
/*----------------------------------------------------------------------------*/

static void
_recv_id_by_src_rank_order(cs_all_to_all_t  *d,
                           const int         src_rank[])
{
  assert(d->recv_id == NULL);

  const cs_lnum_t n_elts = d->n_elts_dest;

  BFT_MALLOC(d->recv_id, n_elts, cs_lnum_t);

  cs_lnum_t n_rs = 0;
  cs_lnum_2_t *rs;

  BFT_MALLOC(rs, n_elts+1, cs_lnum_2_t);

  /* Extract source rank and build rank ranges. */

  int prev_rank = -1;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    if (src_rank[i] != prev_rank) {
      prev_rank = src_rank[i];
      rs[n_rs][0] = src_rank[i];
      rs[n_rs][1] = i;
      n_rs++;
    }
  }

  /* Extend array for future range tests. */

  rs[n_rs][0] = -1;
  rs[n_rs][1] = n_elts;

  /* Order ranges. */

  cs_lnum_t *rs_order;
  BFT_MALLOC(rs_order, n_rs, cs_lnum_t);

  cs_order_lnum_allocated_s(NULL,
                            (const cs_lnum_t  *)rs,
                            2,
                            rs_order,
                            n_rs);

  cs_lnum_t k = 0;
  for (cs_lnum_t i = 0; i < n_rs; i++) {
    cs_lnum_t j = rs_order[i];
    cs_lnum_t s_id = rs[j][1];
    cs_lnum_t e_id = rs[j+1][1];
    for (cs_lnum_t l = s_id; l < e_id; l++) {
      d->recv_id[l] = k;
      k++;
    }
  }
  assert(k == n_elts);

  BFT_FREE(rs_order);
  BFT_FREE(rs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine Crystal Router flags for a distributor.
 *
 * Some flags are only required for the first exchange.
 *
 * \param[in]  d        pointer to associated all-to-all distributor
 * \param[in]  reverse  if true, communicate in reverse direction
 *
 * \return flags required for first Crystal Router call
 */
/*----------------------------------------------------------------------------*/

static int
_cr_flags(cs_all_to_all_t  *d,
          bool              reverse)
{
  cs_assert(d != NULL);

  int cr_flags = 0;

  if (!reverse) {
    if (d->n_elts_dest < 0) {
      if (d->flags & CS_ALL_TO_ALL_USE_DEST_ID)
        cr_flags = cr_flags | CS_CRYSTAL_ROUTER_USE_DEST_ID;
      if (! (d->flags & CS_ALL_TO_ALL_NO_REVERSE)) {
        cr_flags = cr_flags | CS_CRYSTAL_ROUTER_ADD_SRC_ID;
        cr_flags = cr_flags | CS_CRYSTAL_ROUTER_ADD_SRC_RANK;
      }
      if (   d->flags & CS_ALL_TO_ALL_NEED_SRC_RANK
          || d->flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK)
        cr_flags = cr_flags | CS_CRYSTAL_ROUTER_ADD_SRC_RANK;
    }
  }
  else
    cr_flags = cr_flags | CS_CRYSTAL_ROUTER_USE_DEST_ID;

  return cr_flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute destination id by ordering based on source rank after
 *        the first exchange with a Crystal Router.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 * \param[in]  cr  pointer to associated Crystal Router
 */
/*----------------------------------------------------------------------------*/

static void
_cr_recv_id_by_src_rank(cs_all_to_all_t      *d,
                        cs_crystal_router_t  *cr)
{
  cs_assert(d != NULL);

  assert(d->recv_id == NULL);

  int *src_rank;
  BFT_MALLOC(src_rank, d->n_elts_dest_e, int);

  cs_crystal_router_get_data(cr,
                             &src_rank,
                             NULL,
                             NULL,
                             NULL, /* dest_index */
                             NULL);

  if (d->n_elts_dest < 0)
    d->n_elts_dest = cs_crystal_router_n_elts(cr);

  _recv_id_by_src_rank_order(d, src_rank);

  BFT_FREE(src_rank);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if source rank info must be maintained.
 *
 * The source rank is needed if either the CS_CRYSTAL_ROUTER_ADD_SRC_RANK
 * flag is set, or CS_ALL_TO_ALL_NO_REVERSE is not set.
 *
 * It is recommended that if the CS_ALL_TO_ALL_ORDER_BY_SRC_RANK flag is set,
 * for a communication scheme which does not ensure this by default
 * (such as for a Crystal Router), the source rank has been converted
 * to a recv_id, using \ref _recv_id_by_src_rank_order, so if this array
 * is available, the CS_ALL_TO_ALL_ORDER_BY_SRC_RANK flag does imply
 * by itself that source rank info is needed anymore.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return true if source rank info is needed, false otherwise
 */
/*----------------------------------------------------------------------------*/

static bool
_is_src_rank_info_needed(cs_all_to_all_t  *d)
{
  cs_assert(d != NULL);

  bool retval = false;

  if (d->flags & CS_ALL_TO_ALL_NO_REVERSE) {
    if (d->flags & CS_ALL_TO_ALL_NEED_SRC_RANK)
      retval = true;
    else if (d->flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK) {
      if (d->recv_id == NULL && d->n_elts_dest > 0)
        retval = true;
    }
  }
  else
    retval = true;

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor based on destination rank.
 *
 * This is a collective operation on communicator comm.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_USE_DEST_ID,
 * data exchanged will be ordered by the array passed to the
 * \c dest_id argument. For \c n total values received on a rank
 * (as given by \ref cs_all_to_all_n_elts_dest), those destination ids
 * must be in the range [0, \c n[.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
 * data exchanged will be ordered by source rank (this is incompatible
 * with \ref CS_ALL_TO_ALL_USE_DEST_ID.
 *
 * \attention
 * The \c dest_rank and \c dest_id arrays are only referenced by
 * the distributor, not copied, and must remain available throughout
 * the distributor's lifetime. They may be fully transferred to
 * the structure if not needed elsewhere using the
 * \ref cs_all_to_all_transfer_dest_rank and
 * \ref cs_all_to_all_transfer_dest_id functions.
 *
 * \param[in]  n_elts       number of elements
 * \param[in]  flags        sum of ordering and metadata flag constants
 * \param[in]  dest_id      element destination id (required if flags
 *                          contain \ref CS_ALL_TO_ALL_USE_DEST_ID),
 *                          or NULL
 * \param[in]  dest_rank    destination rank for each element
 * \param[in]  comm         associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create(size_t            n_elts,
                     int               flags,
                     const cs_lnum_t  *dest_id,
                     const int         dest_rank[],
                     MPI_Comm          comm)
{
  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  cs_all_to_all_t *d = _all_to_all_create_base(n_elts,
                                               flags,
                                               comm);

  d->dest_id = dest_id;
  d->dest_rank = dest_rank;

  /* Create substructures based on info available at this stage
     (for Crystal Router, delay creation as data is not passed yet) */

  if (d->type == CS_ALL_TO_ALL_MPI_DEFAULT)
    d->dc = _alltoall_caller_create_meta(flags, comm);
  else if (d->type == CS_ALL_TO_ALL_HYBRID)
    d->hc = _hybrid_pex_create_meta(flags, comm);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_TOTAL] += 1;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor for elements whose destination
 * rank is determined from global numbers and block distribution information.
 *
 * This is a collective operation on communicator comm.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_USE_DEST_ID,
 * data exchanged will be ordered by global element number.
 *
 * If the flags bit mask matches \ref CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
 * data exchanged will be ordered by source rank (this is incompatible
 * with \ref CS_ALL_TO_ALL_USE_DEST_ID.
 *
 * \param[in]  n_elts       number of elements
 * \param[in]  flags        sum of ordering and metadata flag constants
 * \param[in]  src_gnum     global source element numbers
 * \param[in]  bi           destination block distribution info
 * \param[in]  comm         associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_from_block(size_t                 n_elts,
                                int                    flags,
                                const cs_gnum_t       *src_gnum,
                                cs_block_dist_info_t   bi,
                                MPI_Comm               comm)
{
  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  cs_all_to_all_t *d = _all_to_all_create_base(n_elts,
                                               flags,
                                               comm);

  /* Compute arrays based on block distribution
     (using global ids and block info at lower levels could
     save some memory, but would make things less modular;
     in addition, rank is an int, while global ids are usually
     64-bit uints, so the overhead is only 50%, and access should
     be faster). */

  BFT_MALLOC(d->_dest_rank, n_elts, int);
  d->dest_rank = d->_dest_rank;

  if (flags & CS_ALL_TO_ALL_USE_DEST_ID) {
    BFT_MALLOC(d->_dest_id, n_elts, cs_lnum_t);
    d->dest_id = d->_dest_id;
  }

  const int rank_step = bi.rank_step;
  const cs_gnum_t block_size = bi.block_size;

  if (d->_dest_id != NULL) {
    #pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (size_t i = 0; i < n_elts; i++) {
      cs_gnum_t g_elt_id = src_gnum[i] -1;
      cs_gnum_t _dest_rank = g_elt_id / block_size;
      d->_dest_rank[i] = _dest_rank*rank_step;
      d->_dest_id[i]   = g_elt_id % block_size;
    }
  }
  else {
    #pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (size_t i = 0; i < n_elts; i++) {
      cs_gnum_t g_elt_id = src_gnum[i] -1;
      cs_gnum_t _dest_rank = g_elt_id / block_size;
      d->_dest_rank[i] = _dest_rank*rank_step;
    }
  }

  /* Create substructures based on info available at this stage
     (for Crystal Router, delay creation as data is not passed yet) */

  if (d->type == CS_ALL_TO_ALL_MPI_DEFAULT)
    d->dc = _alltoall_caller_create_meta(flags, comm);
  else if (d->type == CS_ALL_TO_ALL_HYBRID)
    d->hc = _hybrid_pex_create_meta(flags, comm);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_TOTAL] += 1;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy an all-to-all distributor.
 *
 * \param[in, out]  d   pointer to associated all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_destroy(cs_all_to_all_t **d)
{
  if (d != NULL) {

    cs_timer_t t0, t1;

    t0 = cs_timer_time();

    cs_all_to_all_t *_d = *d;
    if (_d->cr != NULL)
      cs_crystal_router_destroy(&(_d->cr));
    else if (_d->hc != NULL)
      _hybrid_pex_destroy(&(_d->hc));
    else if (_d->dc != NULL)
      _alltoall_caller_destroy(&(_d->dc));

    BFT_FREE(_d->src_rank);
    BFT_FREE(_d->src_id);

    BFT_FREE(_d->_dest_id);
    BFT_FREE(_d->_dest_rank);

    BFT_FREE(_d->recv_id);
    BFT_FREE(_d->src_id);
    BFT_FREE(_d->src_rank);

    BFT_FREE(_d);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                              &t0, &t1);

    /* no increment to _all_to_all_calls[0] as create/destroy are grouped */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer ownership of destination rank to an all-to-all distributor.
 *
 * The dest_rank array should be the same as the one used for the creation of
 * the distributor.
 *
 * \param[in, out]  d          pointer to associated all-to-all distributor
 * \param[in, out]  dest_rank  pointer to element destination rank
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_transfer_dest_rank(cs_all_to_all_t   *d,
                                 int              **dest_rank)
{
  cs_assert(d != NULL);

  if (d->dest_rank == *dest_rank) {
    d->_dest_rank = *dest_rank;
    *dest_rank  = NULL;
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "%s: array transferred (%p)does not match the one used\n"
              "for all-to-all distributor creation (%p).",
              __func__, (void *)*dest_rank, (const void *)d->dest_rank);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer ownership of destination ids to an
 *        all-to-all distributor.
 *
 * The dest_id array should be the same as the one used for the creation of
 * the distributor.
 *
 * \param[in, out]  d        pointer to associated all-to-all distributor
 * \param[in, out]  dest_id  pointer to element destination id
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_transfer_dest_id(cs_all_to_all_t   *d,
                               cs_lnum_t        **dest_id)
{
  cs_assert(d != NULL);

  if (d->dest_id == *dest_id) {
    d->_dest_id = *dest_id;
    *dest_id  = NULL;
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "%s: array transferred (%p)does not match the one used\n"
              "for all-to-all distributor creation (%p).",
              __func__, (void *)*dest_id, (const void *)d->dest_id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements associated with all-to-all distributor.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * If no exchange has been done yet (depending on the communication protocol),
 * metadata will be exchanged by this call, so it is a collective operation.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return  number of elements associated with distributor.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_all_to_all_n_elts_dest(cs_all_to_all_t  *d)
{
  cs_assert(d != NULL);

  /* Obtain count if not available yet */

  if (d->n_elts_dest < 0) {

    cs_timer_t t0, t1;

    t0 = cs_timer_time();

    switch(d->type) {
    case CS_ALL_TO_ALL_MPI_DEFAULT:
      {
        _alltoall_caller_exchange_meta(d->dc,
                                       d->n_elts_src,
                                       d->dest_rank);
        if (d->dc->dest_id_datatype == CS_LNUM_TYPE)
          cs_all_to_all_copy_array(d,
                                   CS_DATATYPE_NULL,
                                   0,
                                   false,
                                   NULL,
                                   NULL);
        else {
          d->n_elts_dest = d->dc->recv_size;
          d->n_elts_dest_e = d->dc->recv_size;
        }
      }
      break;

    case CS_ALL_TO_ALL_HYBRID:
      {
        _hybrid_pex_exchange_meta(d->hc,
                                  d->n_elts_src,
                                  d->dest_rank);
        BFT_FREE(d->_dest_rank);
        if (d->hc->dest_id_datatype == CS_LNUM_TYPE)
          cs_all_to_all_copy_array(d,
                                   CS_DATATYPE_NULL,
                                   0,
                                   false,
                                   NULL,
                                   NULL);
        else {
          d->n_elts_dest = d->hc->recv_size;
          d->n_elts_dest_e = d->hc->recv_size;
        }
      }
      break;

    case CS_ALL_TO_ALL_CRYSTAL_ROUTER:
      {
        cs_crystal_router_t *cr
          = cs_crystal_router_create_s(d->n_elts_src,
                                       0,
                                       CS_DATATYPE_NULL,
                                        _cr_flags(d, false),
                                       NULL,
                                       NULL,
                                       d->dest_id,
                                       d->dest_rank,
                                       d->comm);

        cs_timer_t tcr0 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr0.sec*1e5 + tcr0.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] = 0;
          _all_to_all_trace[_n_trace*9+2] = 0;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_exchange(cr);

        cs_timer_t tcr1 = cs_timer_time();
        cs_timer_counter_add_diff
          (_all_to_all_timers + CS_ALL_TO_ALL_TIME_METADATA, &tcr0, &tcr1);
        _all_to_all_calls[CS_ALL_TO_ALL_TIME_METADATA] += 1;

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr1.sec*1e5 + tcr1.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                            - _all_to_all_trace[(_n_trace-1)*9];
          _all_to_all_trace[_n_trace*9+2] = 1;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        d->n_elts_dest = cs_crystal_router_n_elts(cr);
        d->n_elts_dest_e = cs_crystal_router_n_recv_elts(cr);

        if (d->flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK)
          _cr_recv_id_by_src_rank(d, cr);

        int **p_src_rank = _is_src_rank_info_needed(d) ? &(d->src_rank) : NULL;
        cs_crystal_router_get_data(cr,
                                   p_src_rank,
                                   &(d->recv_id),
                                   &(d->src_id),
                                   NULL, /* dest_index */
                                   NULL);

        cs_crystal_router_destroy(&cr);

      }
      break;

    }

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                              &t0, &t1);
    _all_to_all_calls[CS_ALL_TO_ALL_TIME_TOTAL] += 1;

  }

  return d->n_elts_dest;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate array data using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements returned by \ref cs_all_to_all_n_elts_dest
 * (multiplied by stride and datatype size).
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_data
 * to dest_buffer or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same datatype, stride, and reverse values.
 *
 * \param[in, out]  d          pointer to associated all-to-all distributor
 * \param[in]       datatype   type of data considered
 * \param[in]       stride     number of values per entity (interlaced),
 * \param[in]       reverse    if true, communicate in reverse direction
 * \param[in]       src_data   source data
 * \param[out]      dest_data  pointer to destination data, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

void *
cs_all_to_all_copy_array(cs_all_to_all_t   *d,
                         cs_datatype_t      datatype,
                         int                stride,
                         bool               reverse,
                         const void        *src_data,
                         void              *dest_data)
{
  cs_assert(d != NULL);

  void  *_dest_data = NULL;

  if (_n_trace > 0) {
    fprintf(_all_to_all_trace_bt_log, "\ncs_all_to_all_copy_array: %d\n\n",
            _n_trace);
    cs_base_backtrace_dump(_all_to_all_trace_bt_log, 1);
    bft_printf("cs_all_to_all_copy_array: %d\n", _n_trace);
  }

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Reverse can only be called after direct exchange in most cases
     (this case should be rare, and requires additional echanges,
     but let's play it safe and make sure it works if needed) */

  if (d->n_elts_dest < 0 && reverse)
    cs_all_to_all_copy_array(d,
                             CS_DATATYPE_NULL,
                             0,
                             false,
                             NULL,
                             NULL);

  /* Now do regular exchange */

  switch(d->type) {

  case CS_ALL_TO_ALL_MPI_DEFAULT:
    {
      if (d->n_elts_dest < 0) { /* Exchange metadata if not done yet */
        _alltoall_caller_exchange_meta(d->dc,
                                       d->n_elts_src,
                                       d->dest_rank);
        d->n_elts_dest = d->dc->recv_size;
        d->n_elts_dest_e = d->dc->recv_size;
      }
      size_t n_elts = (reverse) ? d->n_elts_dest : d->n_elts_src;
      if (reverse)
        _alltoall_caller_swap_src_dest(d->dc);
      _alltoall_caller_prepare_s(d->dc,
                                 n_elts,
                                 stride,
                                 datatype,
                                 reverse,
                                 src_data,
                                 d->dest_id,
                                 d->recv_id,
                                 d->dest_rank);
      _dest_data = _alltoall_caller_exchange_s(d,
                                               d->dc,
                                               reverse,
                                               d->dest_rank,
                                               dest_data);
      if (reverse) {
        _alltoall_caller_swap_src_dest(d->dc);
        if (d->dc->send_buffer == src_data)
          d->dc->send_buffer = NULL;
      }
    }
    break;

  case CS_ALL_TO_ALL_HYBRID:
    {
      if (d->n_elts_dest < 0) { /* Exchange metadata if not done yet */
        _hybrid_pex_exchange_meta(d->hc,
                                  d->n_elts_src,
                                  d->dest_rank);
        BFT_FREE(d->_dest_rank);
        d->n_elts_dest = d->hc->recv_size;
        d->n_elts_dest_e = d->hc->recv_size;
      }
      size_t n_elts = (reverse) ? d->n_elts_dest : d->n_elts_src;
      if (reverse)
        _hybrid_pex_swap_src_dest(d->hc);
      _hybrid_pex_prepare_s(d->hc,
                            n_elts,
                            stride,
                            datatype,
                            reverse,
                            src_data,
                            d->dest_id,
                            d->recv_id);
      _dest_data = _hybrid_pex_exchange_s(d,
                                          d->hc,
                                          reverse,
                                          dest_data);
      if (reverse) {
        _hybrid_pex_swap_src_dest(d->hc);
        if (d->hc->send_buffer == src_data)
          d->hc->send_buffer = NULL;
      }
    }
    break;

  case CS_ALL_TO_ALL_CRYSTAL_ROUTER:
    {
      _dest_data = dest_data;
      cs_timer_t tcr0, tcr1;
      cs_crystal_router_t *cr;
      if (!reverse) {
        cr = cs_crystal_router_create_s(d->n_elts_src,
                                        stride,
                                        datatype,
                                        _cr_flags(d, reverse),
                                        src_data,
                                        NULL,
                                        d->dest_id,
                                        d->dest_rank,
                                        d->comm);
        tcr0 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr0.sec*1e5 + tcr0.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] = 0;
          _all_to_all_trace[_n_trace*9+2] = 0;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_exchange(cr);
        tcr1 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          size_t max_sizes[2];
          cs_crystal_router_get_max_sizes(cr, max_sizes);

          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr1.sec*1e5 + tcr1.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                            - _all_to_all_trace[(_n_trace-1)*9];
          _all_to_all_trace[_n_trace*9+2] = 1;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = max_sizes[0];
          _all_to_all_trace[_n_trace*9+6] = max_sizes[1];
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        if (d->n_elts_dest_e < 0) {
          d->n_elts_dest_e = cs_crystal_router_n_recv_elts(cr);
          if (d->flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK)
            _cr_recv_id_by_src_rank(d, cr);
        }
        int **p_src_rank = _is_src_rank_info_needed(d) ? &(d->src_rank) : NULL;
        cs_crystal_router_get_data(cr,
                                   p_src_rank,
                                   &(d->recv_id),
                                   &(d->src_id),
                                   NULL, /* dest_index */
                                   &_dest_data);
        if (d->n_elts_dest < 0)
          d->n_elts_dest = cs_crystal_router_n_elts(cr);
      }
      else {
        const cs_lnum_t *elt_id
          = (d->n_elts_dest < d->n_elts_dest_e) ? d->recv_id : NULL;
        cr = cs_crystal_router_create_s(d->n_elts_dest_e,
                                        stride,
                                        datatype,
                                        _cr_flags(d, reverse),
                                        src_data,
                                        elt_id,
                                        d->src_id,
                                        d->src_rank,
                                        d->comm);
        tcr0 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr0.sec*1e5 + tcr0.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] = 0;
          _all_to_all_trace[_n_trace*9+2] = 0;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_exchange(cr);
        tcr1 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          size_t max_sizes[2];
          cs_crystal_router_get_max_sizes(cr, max_sizes);

          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr1.sec*1e5 + tcr1.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                            - _all_to_all_trace[(_n_trace-1)*9];
          _all_to_all_trace[_n_trace*9+2] = 1;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = max_sizes[0];
          _all_to_all_trace[_n_trace*9+6] = max_sizes[1];
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_get_data(cr,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL, /* dest_index */
                                   &_dest_data);
      }
      cs_crystal_router_destroy(&cr);
      cs_timer_counter_add_diff
        (_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE, &tcr0, &tcr1);
      _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;
    }
    break;

  }

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_TOTAL] += 1;

  return _dest_data;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate local index using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements returned by \ref cs_all_to_all_n_elts_dest.
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_index
 * to dest_index or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same value for the reverse parameter.
 *
 * \param[in, out]  d           pointer to associated all-to-all distributor
 * \param[in]       reverse     if true, communicate in reverse direction
 * \param[in]       src_index   source index
 * \param[out]      dest_index  pointer to destination index, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t *
cs_all_to_all_copy_index(cs_all_to_all_t  *d,
                         bool              reverse,
                         const cs_lnum_t  *src_index,
                         cs_lnum_t        *dest_index)
{
  cs_timer_t t0, t1;

  if (_n_trace > 0) {
    fprintf(_all_to_all_trace_bt_log, "\ncs_all_to_all_copy_index: %d\n\n",
            _n_trace);
    cs_base_backtrace_dump(_all_to_all_trace_bt_log, 1);
    bft_printf("cs_all_to_all_copy_index: %d\n", _n_trace);
  }

  cs_assert(d != NULL);

  cs_lnum_t *src_count = NULL;
  cs_lnum_t *_dest_index = dest_index;

  cs_all_to_all_n_elts_dest(d); /* force sync. if needed */

  cs_lnum_t n_src = (reverse) ? d->n_elts_dest : d->n_elts_src;
  cs_lnum_t n_dest = -1;

  if (dest_index == NULL) {
    n_dest = (reverse) ? d->n_elts_src : d->n_elts_dest;
  }

  t0 = cs_timer_time();

  if (dest_index == NULL)
    BFT_MALLOC(_dest_index, n_dest + 1, cs_lnum_t);

  /* Convert send index to count, then exchange */

  BFT_MALLOC(src_count, n_src, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_src; i++)
    src_count[i] = src_index[i+1] - src_index[i];

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);

  cs_all_to_all_copy_array(d,
                           CS_LNUM_TYPE,
                           1,
                           reverse,
                           src_count,
                           _dest_index + 1);

  t0 = cs_timer_time();

  BFT_FREE(src_count);

  _dest_index[0] = 0;

  if (n_dest < 1)
    n_dest = (reverse) ? d->n_elts_src : d->n_elts_dest;

  for (cs_lnum_t i = 0; i < n_dest; i++)
    _dest_index[i+1] += _dest_index[i];

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);

  return _dest_index;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Communicate local index using all-to-all distributor.
 *
 * If a destination buffer is provided, it should be of sufficient size for
 * the number of elements indicated by
 * dest_index[\ref cs_all_to_all_n_elts_dest "cs_all_to_all_n_elts_dest(d)"];
 *
 * If no buffer is provided, one is allocated automatically, and transferred
 * to the caller (who is responsible for freeing it when no longer needed).
 *
 * If used in reverse mode, data is still communicated from src_index
 * to dest_index or an internal buffer, but communication direction
 * (i.e. source and destination ranks) are reversed.
 *
 * This is obviously a collective operation, and all ranks must provide
 * the same value for the reverse parameter.
 *
 * \param[in, out]  d           pointer to associated all-to-all distributor
 * \param[in]       datatype    type of data considered
 * \param[in]       reverse     if true, communicate in reverse direction
 * \param[in]       src_index   source index
 * \param[in]       src_data    source data
 * \param[in]       dest_index  destination index
 * \param[out]      dest_data   pointer to destination data, or NULL
 *
 * \return pointer to destination data (dest_buffer if non-NULL)
 */
/*----------------------------------------------------------------------------*/

void *
cs_all_to_all_copy_indexed(cs_all_to_all_t  *d,
                           cs_datatype_t     datatype,
                           bool              reverse,
                           const cs_lnum_t  *src_index,
                           const void       *src_data,
                           const cs_lnum_t  *dest_index,
                           void             *dest_data)
{
  cs_assert(d != NULL);

  if (_n_trace > 0) {
    fprintf(_all_to_all_trace_bt_log, "\ncs_all_to_all_copy_indexed: %d\n\n",
            _n_trace);
    cs_base_backtrace_dump(_all_to_all_trace_bt_log, 1);
    bft_printf("cs_all_to_all_copy_indexed: %d\n", _n_trace);
  }

  void  *_dest_data = NULL;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Reverse can only be called after direct exchange in most cases
     (this case should be rare, and requires additional echanges,
     but let's play it safe and make sure it works if needed) */

  if (d->n_elts_dest < 0 && reverse)
    cs_all_to_all_copy_array(d,
                             CS_DATATYPE_NULL,
                             0,
                             false,
                             NULL,
                             NULL);

  /* Now do regular exchange */

  switch(d->type) {

  case CS_ALL_TO_ALL_MPI_DEFAULT:
    {
      if (d->n_elts_dest < 0) { /* Exchange metadata if not done yet */
        _alltoall_caller_exchange_meta(d->dc,
                                       d->n_elts_src,
                                       d->dest_rank);
        d->n_elts_dest = d->dc->recv_size;
      }
      size_t n_elts = (reverse) ? d->n_elts_dest : d->n_elts_src;
      _alltoall_caller_save_meta_i(d->dc);
      if (reverse)
        _alltoall_caller_swap_src_dest(d->dc);
      _alltoall_caller_prepare_i(d->dc,
                                 n_elts,
                                 datatype,
                                 reverse,
                                 src_index,
                                 dest_index,
                                 src_data,
                                 d->recv_id,
                                 d->dest_rank);
      _dest_data = _alltoall_caller_exchange_i(d,
                                               d->dc,
                                               reverse,
                                               dest_index,
                                               d->dest_rank,
                                               dest_data);
      if (reverse) {
        _alltoall_caller_swap_src_dest(d->dc);
        if (d->dc->send_buffer == src_data)
          d->dc->send_buffer = NULL;
      }
      _alltoall_caller_reset_meta_i(d->dc, d->dest_rank);
    }
    break;

  case CS_ALL_TO_ALL_HYBRID:
    {
      if (d->n_elts_dest < 0) { /* Exchange metadata if not done yet */
        _hybrid_pex_exchange_meta(d->hc,
                                  d->n_elts_src,
                                  d->dest_rank);
        BFT_FREE(d->_dest_rank);
        d->n_elts_dest = d->hc->recv_size;
      }
      size_t n_elts = (reverse) ? d->n_elts_dest : d->n_elts_src;
      _hybrid_pex_save_meta_i(d->hc);
      if (reverse)
        _hybrid_pex_swap_src_dest(d->hc);
      _hybrid_pex_prepare_i(d->hc,
                            n_elts,
                            datatype,
                            reverse,
                            src_index,
                            dest_index,
                            src_data,
                            d->recv_id);
      _dest_data = _hybrid_pex_exchange_i(d,
                                          d->hc,
                                          reverse,
                                          dest_index,
                                          dest_data);
      if (reverse) {
        _hybrid_pex_swap_src_dest(d->hc);
        if (d->hc->send_buffer == src_data)
          d->hc->send_buffer = NULL;
      }
      _hybrid_pex_reset_meta_i(d->hc);
    }
    break;

  case CS_ALL_TO_ALL_CRYSTAL_ROUTER:
    {
      _dest_data = dest_data;
      cs_timer_t tcr0, tcr1;
      cs_crystal_router_t *cr;
      if (!reverse) {
        cr = cs_crystal_router_create_i(d->n_elts_src,
                                        datatype,
                                        _cr_flags(d, reverse),
                                        src_index,
                                        src_data,
                                        NULL,
                                        d->dest_id,
                                        d->dest_rank,
                                        d->comm);
        tcr0 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr0.sec*1e5 + tcr0.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] = 0;
          _all_to_all_trace[_n_trace*9+2] = 0;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_exchange(cr);
        tcr1 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          size_t max_sizes[2];
          cs_crystal_router_get_max_sizes(cr, max_sizes);

          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr1.sec*1e5 + tcr1.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                            - _all_to_all_trace[(_n_trace-1)*9];
          _all_to_all_trace[_n_trace*9+2] = 1;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = max_sizes[0];
          _all_to_all_trace[_n_trace*9+6] = max_sizes[1];
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        if (d->n_elts_dest < 0) {
          d->n_elts_dest = cs_crystal_router_n_elts(cr);
          d->n_elts_dest_e = cs_crystal_router_n_recv_elts(cr);
        }
        int **p_src_rank = (d->src_rank == NULL) ? &(d->src_rank) : NULL;
        cs_crystal_router_get_data(cr,
                                   p_src_rank,
                                   &(d->recv_id),
                                   &(d->src_id),
                                   NULL, /* dest_index */
                                   &_dest_data);
      }
      else {
        const cs_lnum_t *elt_id
          = (d->n_elts_dest < d->n_elts_dest_e) ? d->recv_id : NULL;
        cr = cs_crystal_router_create_i(d->n_elts_dest_e,
                                        datatype,
                                        _cr_flags(d, reverse),
                                        src_index,
                                        src_data,
                                        elt_id,
                                        d->src_id,
                                        d->src_rank,
                                        d->comm);
        tcr0 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr0.sec*1e5 + tcr0.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] = 0;
          _all_to_all_trace[_n_trace*9+2] = 0;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = 0;
          _all_to_all_trace[_n_trace*9+6] = 0;
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_exchange(cr);
        tcr1 = cs_timer_time();

        if (_n_trace < _n_trace_max) {
          size_t max_sizes[2];
          cs_crystal_router_get_max_sizes(cr, max_sizes);

          /* Time to 1-5 s */
          _all_to_all_trace[_n_trace*9] = tcr1.sec*1e5 + tcr1.nsec/1e4;
          _all_to_all_trace[_n_trace*9+1] =   _all_to_all_trace[_n_trace*9]
                                            - _all_to_all_trace[(_n_trace-1)*9];
          _all_to_all_trace[_n_trace*9+2] = 1;
          _all_to_all_trace[_n_trace*9+3] = bft_mem_usage_pr_size();
          _all_to_all_trace[_n_trace*9+4] = bft_mem_usage_max_pr_size();
          _all_to_all_trace[_n_trace*9+5] = max_sizes[0];
          _all_to_all_trace[_n_trace*9+6] = max_sizes[1];
          _all_to_all_trace[_n_trace*9+7] = 0;
          _all_to_all_trace[_n_trace*9+8] = 0;
          _n_trace += 1;
        }

        cs_crystal_router_get_data(cr,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL, /* dest_index */
                                   &_dest_data);
      }
      cs_crystal_router_destroy(&cr);
      cs_timer_counter_add_diff
        (_all_to_all_timers + CS_ALL_TO_ALL_TIME_EXCHANGE, &tcr0, &tcr1);
      _all_to_all_calls[CS_ALL_TO_ALL_TIME_EXCHANGE] += 1;
    }
    break;

  }

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);
  _all_to_all_calls[CS_ALL_TO_ALL_TIME_TOTAL] += 1;

  return _dest_data;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get array of source element ranks associated with an
 *        all-to-all distributor.
 *
 * This function should be called only after \ref cs_all_to_all_copy_array,
 * and allocates and returns an array of source element ranks matching the
 * exchanged data elements.
 *
 * It should also be called only if the distributor creation flags match
 * CS_ALL_TO_ALL_NEED_SRC_RANK or CS_ALL_TO_ALL_ORDER_BY_SRC_RANK.
 *
 * The returned data is owned by the caller, who is responsible for freeing
 * it when no longer needed.
 *
 * If source ranks are not available (depending on the distributor
 * creation function and options), the matching pointer will
 * be set to NULL.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return  array of source ranks (or NULL)
 */
/*----------------------------------------------------------------------------*/

int *
cs_all_to_all_get_src_rank(cs_all_to_all_t  *d)
{
  cs_timer_t t0 = cs_timer_time();

  int *src_rank = NULL;

  cs_assert(d != NULL);

  if (! (   d->flags & CS_ALL_TO_ALL_NEED_SRC_RANK
         || d->flags & CS_ALL_TO_ALL_ORDER_BY_SRC_RANK))
    bft_error(__FILE__, __LINE__, 0,
              "%s: is called for a distributor with flags %d, which does not\n"
              "match masks CS_ALL_TO_ALL_NEED_SRC_RANK (%d) or "
              "CS_ALL_TO_ALL_ORDER_BY_SRC_RANK (%d).",
              __func__, d->flags,
              CS_ALL_TO_ALL_NEED_SRC_RANK,
              CS_ALL_TO_ALL_ORDER_BY_SRC_RANK);

  BFT_MALLOC(src_rank, d->n_elts_dest, int);

  switch(d->type) {

  case CS_ALL_TO_ALL_MPI_DEFAULT:
    {
      const _mpi_all_to_all_caller_t *dc = d->dc;

      for (int i = 0; i < dc->n_ranks; i++) {
        for (cs_lnum_t j = dc->recv_displ[i]; j < dc->recv_displ[i+1]; j++)
          src_rank[j] = i;
      }
    }
    break;

  case CS_ALL_TO_ALL_HYBRID:
    {
      const _hybrid_pex_t *hc = d->hc;

      for (int i = 0; i < hc->rn_recv->size; i++) {
        int rank_id = hc->rn_recv->rank[i];
        for (cs_lnum_t j = hc->recv_displ[i]; j < hc->recv_displ[i+1]; j++)
          src_rank[j] = rank_id;
      }
    }
    break;

  case CS_ALL_TO_ALL_CRYSTAL_ROUTER:
    {
      if (d->src_rank != NULL)
        memcpy(src_rank, d->src_rank, d->n_elts_dest*sizeof(int));
    }

  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers + CS_ALL_TO_ALL_TIME_TOTAL,
                            &t0, &t1);

  return src_rank;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of all-to-all distributor algorithm choice.
 *
 * \return  current type of all-to-all distributor algorithm choice
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_type_t
cs_all_to_all_get_type(void)
{
  return _all_to_all_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current type of hybrid all-to-all distributor parameters.
 *
 * \param[out]  rne_type  type of metadata exchange algorithm, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_hybrid_parameters(cs_rank_neighbors_exchange_t  *rne_type)
{
  if (rne_type != NULL)
    *rne_type = _hybrid_meta_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of all-to-all distributor algorithm choice.
 *
 * \param[in]  rne_type  type of metadata exchange algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_set_hybrid_parameters(cs_rank_neighbors_exchange_t  rne_type)
{
  _hybrid_meta_type = rne_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set current type of all-to-all distributor algorithm choice.
 *
 * \param  t  type of all-to-all distributor algorithm choice to select
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_set_type(cs_all_to_all_type_t  t)
{
  _all_to_all_type = t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log performance information relative to instrumented all-to-all
 * distribution.
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_log_finalize(void)
{
  cs_crystal_router_log_finalize();

#if defined(HAVE_MPI)

  if (_all_to_all_calls[0] <= 0)
    return;

  char method_name[96];
  switch(_all_to_all_type) {
  case CS_ALL_TO_ALL_MPI_DEFAULT:
    snprintf(method_name, 96, N_("MPI_Alltoall and MPI_Alltoallv"));
    break;
  case CS_ALL_TO_ALL_HYBRID:
    snprintf(method_name, 96, N_("Hybrid, %s (metadata), %s (data)"),
             _(cs_rank_neighbors_exchange_name[_hybrid_meta_type]),
             "MPI_Alltoallv");
    break;
  case CS_ALL_TO_ALL_CRYSTAL_ROUTER:
    snprintf(method_name, 96, N_("Crystal Router algorithm"));
    break;
  }
  method_name[95] = '\0';

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nAll-to-many operations: (%s):\n\n"),
                method_name);

  /* Print times */

  double wtimes[3] = {0, 0, 0};
  double wtimes_mean[3], wtimes_max[3], wtimes_min[3];
  for (int i = 0; i < 3; i++) {
    if (_all_to_all_calls[i] > 0)
      wtimes[i] = (_all_to_all_timers[i]).nsec*1e-9;
    wtimes_mean[i] = wtimes[i];
    wtimes_max[i] = wtimes[i];
    wtimes_min[i] = wtimes[i];
  }
  if (cs_glob_n_ranks > 1) {
    MPI_Allreduce(wtimes, wtimes_mean, 3, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(wtimes, wtimes_max, 3, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
    MPI_Allreduce(wtimes, wtimes_min, 3, MPI_DOUBLE, MPI_MIN, cs_glob_mpi_comm);
  }
  for (int i = 0; i < 3; i++)
    wtimes_mean[i] /= cs_glob_n_ranks;

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("                             mean        minimum      maximum"
       "     calls\n"
       "  Total:             %12.5f s %12.5f %12.5f s   %lu\n"
       "  Metadata exchange: %12.5f s %12.5f %12.5f s   %lu\n"
       "  Data exchange:     %12.5f s %12.5f %12.5f s   %lu\n\n"),
     wtimes_mean[0], wtimes_min[0], wtimes_max[0],
     (unsigned long)(_all_to_all_calls[0]),
     wtimes_mean[1], wtimes_min[1], wtimes_max[1],
     (unsigned long)(_all_to_all_calls[1]),
     wtimes_mean[2], wtimes_min[2], wtimes_max[2],
     (unsigned long)(_all_to_all_calls[2]));

  cs_log_separator(CS_LOG_PERFORMANCE);

  if (cs_glob_n_ranks > 1 && _n_trace > 0) {
    cs_gnum_t *_all_to_all_sum;
    BFT_MALLOC(_all_to_all_sum, _n_trace*9, uint64_t);
    cs_gnum_t *_all_to_all_max;
    BFT_MALLOC(_all_to_all_max, _n_trace*9, uint64_t);

    MPI_Allreduce(_all_to_all_trace, _all_to_all_sum, _n_trace*9, MPI_UINT64_T,
                  MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(_all_to_all_trace, _all_to_all_max, _n_trace*9, MPI_UINT64_T,
                  MPI_MAX, cs_glob_mpi_comm);

    for (int j = 0; j < _n_trace*9; j++)
      _all_to_all_sum[j] /= cs_glob_n_ranks;

    if (cs_glob_rank_id < 1) {
      FILE *f = fopen("all_to_all_trace.csv", "w");
      fprintf(f, "call, time, dt_mean, dt_max, stage, mem_cur_mean, mem_cur_max, "
              "mem_max_mean, mem_max, send_size_mean, send_size_max, "
              "recv_size_mean, recv_size_max, "
              "send_ranks_mean, send_ranks_max, "
              "recv_ranks_mean, recv_ranks_max\n");
      for (int j = 0; j < _n_trace; j++)
        fprintf(f,
                "%d, %g, %g, %g, %d,"
                "%llu, %llu, %llu, %llu, %llu, %llu,"
                "%llu, %llu, %llu, %llu, %llu, %llu\n",
                j,
                (_all_to_all_trace[j*9] - _all_to_all_trace[0])/1.e5,
                _all_to_all_sum[j*9+1]/1.e5,
                _all_to_all_max[j*9+1]/1.e5,
                (int)_all_to_all_trace[j*9+2],
                (unsigned long long)_all_to_all_sum[j*9+3],
                (unsigned long long)_all_to_all_max[j*9+3],
                (unsigned long long)_all_to_all_sum[j*9+4],
                (unsigned long long)_all_to_all_max[j*9+4],
                (unsigned long long)_all_to_all_sum[j*9+5]/1000,
                (unsigned long long)_all_to_all_max[j*9+5]/1000,
                (unsigned long long)_all_to_all_sum[j*9+6]/1000,
                (unsigned long long)_all_to_all_max[j*9+6]/1000,
                (unsigned long long)_all_to_all_sum[j*9+7],
                (unsigned long long)_all_to_all_max[j*9+7],
                (unsigned long long)_all_to_all_sum[j*9+8],
                (unsigned long long)_all_to_all_max[j*9+8]);
      fclose(f);
      fclose(_all_to_all_trace_bt_log);
      _all_to_all_trace_bt_log = NULL;
    }

    BFT_FREE(_all_to_all_sum);
    BFT_FREE(_all_to_all_max);
    BFT_FREE(_all_to_all_trace);
  }


#endif /* defined(HAVE_MPI) */
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
