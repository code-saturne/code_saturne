/*============================================================================
 * All-to-all parallel data exchange.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_block_dist.h"
#include "cs_log.h"
#include "cs_order.h"
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
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Metadata flags
 */

#define CS_ALL_TO_ALL_USE_GNUM           (1 << 0)
#define CS_ALL_TO_ALL_USE_DEST_ID        (1 << 1)
#define CS_ALL_TO_ALL_USE_SRC_ID         (1 << 2)

/*=============================================================================
 * Local type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

typedef struct {

  cs_datatype_t   datatype;          /* associated datatype */
  cs_datatype_t   dest_id_datatype;  /* type of destination id (CS_GNUM_TYPE,
                                        CS_LNUM_TYPE or CS_DATATYPE_NULL) */
  bool            add_src_id;        /* add source id ? */

  size_t          stride;            /* stride if strided, 0 otherwise */

  size_t          src_id_shift;      /* starting byte for source id */
  size_t          elt_shift;         /* starting byte for element data */
  size_t          comp_size;         /* Composite element size, with padding */

  size_t          send_size;
  size_t          recv_size;

  unsigned char  *send_buffer;
  unsigned char  *recv_buffer;

  int            *send_count;        /* Send counts for MPI_Alltoall */
  int            *recv_count;        /* Receive counts for MPI_Alltoall */
  int            *send_displ;        /* Send displs for MPI_Alltoall */
  int            *recv_displ;        /* Receive displs for MPI_Alltoall */

  int            *src_rank;          /* Source rank, or NULL */
  int            *dest_rank;         /* Destination rank, or NULL */

  MPI_Comm        comm;              /* Associated MPI communicator */
  MPI_Datatype    comp_type;         /* Associated MPI datatype */
  int             rank_id;           /* local rank id in comm */

  int             n_ranks;           /* Number of ranks associated with
                                        communicator */

} _mpi_all_to_all_caller_t;

/* Crystal router management structure */

typedef struct { /* Structure used to manage crystal router information */

  cs_datatype_t   datatype;          /* associated datatype */
  cs_datatype_t   dest_id_datatype;  /* type of destination id (CS_GNUM_TYPE,
                                        CS_LNUM_TYPE or CS_DATATYPE_NULL) */
  bool            add_src_id;        /* add source id ? */

  size_t          stride;            /* stride if strided, 0 otherwise */

  size_t          dest_id_shift;     /* starting byte for destination id */
  size_t          src_id_shift;      /* starting byte for source id */
  size_t          elt_shift;         /* starting byte for element data */

  size_t          elt_size;          /* element size if strided, 0 otherwise */
  size_t          comp_size;         /* composite metadata + element size if
                                        strided, 0 otherwise */
  size_t          n_elts[2];
  size_t          buffer_size[2];
  unsigned char  *buffer[2];

  MPI_Comm        comm;              /* associated MPI communicator */
  MPI_Datatype    comp_type;         /* Associated MPI datatype */
  int             rank_id;           /* local rank id in comm */
  int             n_ranks;           /* comm size */

} _crystal_router_t;

#endif /* defined(HAVE_MPI) */

/* Structure used to redistribute data */

#if defined(HAVE_MPI)

struct _cs_all_to_all_t {

  bool                       strided;  /* True if strided, false otherwise */

  _mpi_all_to_all_caller_t  *dc;       /* Default MPI_Alltoall(v) caller */
  _crystal_router_t         *cr;       /* associated crystal-router */

};

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default all-to-all type */

static cs_all_to_all_type_t _all_to_all_type = CS_ALL_TO_ALL_MPI_DEFAULT;

#if defined(HAVE_MPI)

/* Call counter and timer: 0: setup, 1: exchange,
   2: swap source and destination, 3: sort by source rank, 4: copy data */

static size_t              _all_to_all_calls[5] = {0, 0, 0, 0, 0};
static cs_timer_counter_t  _all_to_all_timers[5];

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(HAVE_MPI)

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

  for (i = 1; i < n_ranks; i++)
    displ[i] = displ[i-1] + count[i-1];

  total_count = displ[n_ranks-1] + count[n_ranks-1];

  return total_count;
}

/*----------------------------------------------------------------------------
 * First stage of creation for an MPI_Alltoall(v) caller for strided data.
 *
 * parameters:
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL depending on elt_id values)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _mpi_all_to_all_caller_t *
_alltoall_caller_create_meta_s(int            stride,
                               cs_datatype_t  datatype,
                               cs_datatype_t  dest_id_datatype,
                               bool           add_src_id,
                               MPI_Comm       comm)
{
  int rank_id, n_ranks;
  size_t elt_size = cs_datatype_size[datatype]*stride;
  size_t align_size = sizeof(cs_lnum_t);

  _mpi_all_to_all_caller_t *dc = NULL;

  /* Communicator info */

  MPI_Comm_rank(comm, &rank_id);
  MPI_Comm_size(comm, &n_ranks);

  /* Allocate structure */

  BFT_MALLOC(dc, 1, _mpi_all_to_all_caller_t);

  if (datatype == CS_GNUM_TYPE || datatype == CS_LNUM_TYPE)
    dc->datatype = datatype;
  else
    dc->datatype = CS_DATATYPE_NULL;

  dc->dest_id_datatype = dest_id_datatype;
  dc->add_src_id = add_src_id;

  dc->stride = stride;

  dc->send_size = 0;
  dc->recv_size = 0;

  dc->send_size = 0;
  dc->recv_size = 0;

  dc->src_rank = NULL;
  dc->dest_rank = NULL;

  dc->comm = comm;
  dc->rank_id = rank_id;
  dc->n_ranks = n_ranks;

  dc->send_buffer = NULL;
  dc->recv_buffer = NULL;

  BFT_MALLOC(dc->send_count, dc->n_ranks, int);
  BFT_MALLOC(dc->recv_count, dc->n_ranks, int);
  BFT_MALLOC(dc->send_displ, dc->n_ranks + 1, int);
  BFT_MALLOC(dc->recv_displ, dc->n_ranks + 1, int);

  /* Compute data size and alignment */

  if (dc->dest_id_datatype == CS_GNUM_TYPE) {
    dc->src_id_shift = sizeof(cs_gnum_t);
    align_size = sizeof(cs_gnum_t);
  }
  else if (dc->dest_id_datatype == CS_LNUM_TYPE)
    dc->src_id_shift = sizeof(cs_lnum_t);
  else
    dc->src_id_shift = 0;

  if (dc->add_src_id)
    dc->elt_shift = dc->src_id_shift + cs_datatype_size[CS_LNUM_TYPE];
  else
    dc->elt_shift = dc->src_id_shift;

  if (dc->elt_shift % align_size)
    dc->elt_shift += align_size - (dc->elt_shift % align_size);

  dc->comp_size = dc->elt_shift + elt_size;

  if (elt_size % align_size)
    dc->comp_size += align_size - (elt_size % align_size);;

  /* Create associated MPI datatype */

  MPI_Type_contiguous(dc->comp_size, MPI_BYTE, &(dc->comp_type));
  MPI_Type_commit(&(dc->comp_type));

  /* Return pointer to structure */

  return dc;
}

/*----------------------------------------------------------------------------
 * Create a MPI_Alltoall(v) caller for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL depending on elt_id values)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   dest_id          <-- element destination id, global number, or NULL
 *   dest_rank        <-- destination rank for each element
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _mpi_all_to_all_caller_t *
_alltoall_caller_create_s(size_t         n_elts,
                          int            stride,
                          cs_datatype_t  datatype,
                          cs_datatype_t  dest_id_datatype,
                          bool           add_src_id,
                          void          *elt,
                          void          *elt_id,
                          const int      dest_rank[],
                          MPI_Comm       comm)
{
  int i;
  size_t j, k;

  size_t elt_size = cs_datatype_size[datatype]*stride;

  unsigned char *_elt = elt;
  unsigned char *_elt_id = elt_id;

  assert(elt != NULL || n_elts == 0);

  /* Create base data */

  _mpi_all_to_all_caller_t *dc = _alltoall_caller_create_meta_s(stride,
                                                                datatype,
                                                                dest_id_datatype,
                                                                add_src_id,
                                                                comm);

  /* Count values to send and receive */

  for (i = 0; i < dc->n_ranks; i++)
    dc->send_count[i] = 0;

  for (j = 0; j < n_elts; j++)
    dc->send_count[dest_rank[j]] += 1;

  dc->send_size = _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);

  /* Allocate send buffer */

  BFT_MALLOC(dc->send_buffer, dc->send_size*dc->comp_size, unsigned char);
  memset(dc->send_buffer, 0, dc->send_size*dc->comp_size);

  /* Copy metadata */

  if (dc->dest_id_datatype != CS_DATATYPE_NULL) {
    const size_t id_size =   (dc->dest_id_datatype == CS_GNUM_TYPE)
                           ? sizeof(cs_gnum_t) : sizeof(cs_lnum_t);
    for (j = 0; j < n_elts; j++) {
      size_t w_displ = dc->send_displ[dest_rank[j]]*dc->comp_size;
      size_t r_displ = j*id_size;
      dc->send_displ[dest_rank[j]] += 1;
      for (k = 0; k < id_size; k++)
        dc->send_buffer[w_displ + k] = _elt_id[r_displ + k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }

  if (dc->add_src_id) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (j = 0; j < n_elts; j++) {
      cs_lnum_t src_id = j;
      unsigned char *_src_id = (unsigned char *)(&src_id);
      size_t w_displ =   dc->send_displ[dest_rank[j]]*dc->comp_size
                       + dc->src_id_shift;
      dc->send_displ[dest_rank[j]] += 1;
      for (k = 0; k < id_size; k++)
        dc->send_buffer[w_displ + k] = _src_id[k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }

  /* Copy data */

  for (j = 0; j < n_elts; j++) {

    size_t w_displ =   dc->send_displ[dest_rank[j]]*dc->comp_size
                     + dc->elt_shift;
    size_t r_displ = j*elt_size;

    dc->send_displ[dest_rank[j]] += 1;

    for (k = 0; k < elt_size; k++)
      dc->send_buffer[w_displ + k] = _elt[r_displ + k];

  }
  /* Reset send_displ */
  for (i = 0; i < dc->n_ranks; i++)
    dc->send_displ[i] -= dc->send_count[i];

  /* Return pointer to structure */

  return dc;
}

/*----------------------------------------------------------------------------
 * Create a MPI_Alltoall(v) caller for strided data using block information.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   elt_gnum         <-- global element numbers
 *   bi               <-- destination block distribution info
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _mpi_all_to_all_caller_t *
_alltoall_caller_create_from_block_s(size_t                 n_elts,
                                     int                    stride,
                                     cs_datatype_t          datatype,
                                     cs_datatype_t          dest_id_datatype,
                                     bool                   add_src_id,
                                     void                  *elt,
                                     const cs_gnum_t       *elt_gnum,
                                     cs_block_dist_info_t   bi,
                                     MPI_Comm               comm)
{
  int i;
  size_t j, k;

  size_t elt_size = cs_datatype_size[datatype]*stride;

  unsigned char *_elt = elt;

  const int rank_step = bi.rank_step;
  const cs_gnum_t block_size = bi.block_size;

  assert(elt != NULL || n_elts == 0);

  /* Create base data */

  _mpi_all_to_all_caller_t *dc = _alltoall_caller_create_meta_s(stride,
                                                                datatype,
                                                                dest_id_datatype,
                                                                add_src_id,
                                                                comm);

  /* Count values to send and receive */

  for (i = 0; i < dc->n_ranks; i++)
    dc->send_count[i] = 0;

  for (j = 0; j < n_elts; j++) {
    cs_gnum_t g_elt_id = elt_gnum[j] -1;
    cs_gnum_t _dest_rank = g_elt_id / block_size;
    int dest_rank = _dest_rank*rank_step;
    dc->send_count[dest_rank] += 1;
  }

  dc->send_size = _compute_displ(dc->n_ranks, dc->send_count, dc->send_displ);

  /* Allocate send buffer */

  BFT_MALLOC(dc->send_buffer, dc->send_size*dc->comp_size, unsigned char);
  memset(dc->send_buffer, 0, dc->send_size*dc->comp_size);

  /* Copy metadata */

  if (dc->dest_id_datatype == CS_GNUM_TYPE) {
    const size_t id_size = sizeof(cs_gnum_t);
    for (j = 0; j < n_elts; j++) {
      cs_gnum_t g_elt_num = elt_gnum[j];
      cs_gnum_t g_elt_id = g_elt_num -1;
      cs_gnum_t _dest_rank = g_elt_id / block_size;
      int dest_rank = _dest_rank*rank_step;
      size_t w_displ = dc->send_displ[dest_rank]*dc->comp_size;
      unsigned char *_g_elt_num = (unsigned char *)(&g_elt_num);
      dc->send_displ[dest_rank] += 1;
      for (k = 0; k < id_size; k++)
        dc->send_buffer[w_displ + k] = _g_elt_num[k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }
  else if (dc->dest_id_datatype == CS_LNUM_TYPE) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (j = 0; j < n_elts; j++) {
      cs_gnum_t g_elt_id = elt_gnum[j] -1;
      cs_gnum_t _dest_rank = g_elt_id / block_size;
      cs_lnum_t b_elt_id = g_elt_id % block_size;;
      int dest_rank = _dest_rank*rank_step;
      size_t w_displ = dc->send_displ[dest_rank]*dc->comp_size;
      unsigned char *_b_elt_id = (unsigned char *)(&b_elt_id);
      dc->send_displ[dest_rank] += 1;
      for (k = 0; k < id_size; k++)
        dc->send_buffer[w_displ + k] = _b_elt_id[k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }

  if (dc->add_src_id) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (j = 0; j < n_elts; j++) {
      cs_gnum_t g_elt_id = elt_gnum[j] -1;
      cs_gnum_t _dest_rank = g_elt_id / block_size;
      int dest_rank = _dest_rank*rank_step;
      cs_lnum_t src_id = j;
      unsigned char *_src_id = (unsigned char *)(&src_id);
      size_t w_displ =   dc->send_displ[dest_rank]*dc->comp_size
                       + dc->src_id_shift;
      dc->send_displ[dest_rank] += 1;
      for (k = 0; k < id_size; k++)
        dc->send_buffer[w_displ + k] = _src_id[k];
    }
    /* Reset send_displ */
    for (i = 0; i < dc->n_ranks; i++)
      dc->send_displ[i] -= dc->send_count[i];
  }

  /* Copy data */

  for (j = 0; j < n_elts; j++) {

    cs_gnum_t g_elt_id = elt_gnum[j] -1;
    cs_gnum_t _dest_rank = g_elt_id / block_size;
    int dest_rank = _dest_rank*rank_step;

    size_t w_displ =   dc->send_displ[dest_rank]*dc->comp_size
                     + dc->elt_shift;
    size_t r_displ = j*elt_size;

    dc->send_displ[dest_rank] += 1;

    for (k = 0; k < elt_size; k++)
      dc->send_buffer[w_displ + k] = _elt[r_displ + k];

  }
  /* Reset send_displ */
  for (i = 0; i < dc->n_ranks; i++)
    dc->send_displ[i] -= dc->send_count[i];

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
    BFT_FREE(_dc->recv_buffer);
    BFT_FREE(_dc->send_buffer);
    BFT_FREE(_dc->dest_rank);
    BFT_FREE(_dc->src_rank);
    BFT_FREE(_dc->recv_displ);
    BFT_FREE(_dc->send_displ);
    BFT_FREE(_dc->recv_count);
    BFT_FREE(_dc->send_count);
    BFT_FREE(*dc);
  }
}

/*----------------------------------------------------------------------------
 * Exchange strided data with a MPI_Alltoall(v) caller.
 *
 * Order of data from a same source rank is preserved.
 *
 * parameters:
 *   dc <-> associated MPI_Alltoall(v) caller structure
 *---------------------------------------------------------------------------*/

static  void
_alltoall_caller_exchange_strided(_mpi_all_to_all_caller_t  *dc)
{
  /* Exchange counts */

  MPI_Alltoall(dc->send_count, 1, MPI_INT,
               dc->recv_count, 1, MPI_INT,
               dc->comm);

  dc->recv_size = _compute_displ(dc->n_ranks, dc->recv_count, dc->recv_displ);

  /* Exchange data */

  BFT_REALLOC(dc->recv_buffer, dc->recv_size*dc->comp_size, unsigned char);

  MPI_Alltoallv(dc->send_buffer, dc->send_count, dc->send_displ, dc->comp_type,
                dc->recv_buffer, dc->recv_count, dc->recv_displ, dc->comp_type,
                dc->comm);
}

/*----------------------------------------------------------------------------
 * First stage of creation for a crystal router for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL depending on elt_id values)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _crystal_router_t *
_crystal_create_meta_s(size_t         n_elts,
                       int            stride,
                       cs_datatype_t  datatype,
                       cs_datatype_t  dest_id_datatype,
                       bool           add_src_id,
                       MPI_Comm       comm)
{
  int rank_id, n_ranks;
  _crystal_router_t *cr = NULL;

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

  BFT_MALLOC(cr, 1, _crystal_router_t);

  if (datatype == CS_GNUM_TYPE || datatype == CS_LNUM_TYPE)
    cr->datatype = datatype;
  else
    cr->datatype = CS_DATATYPE_NULL;

  cr->dest_id_datatype = dest_id_datatype;
  cr->add_src_id = add_src_id;

  cr->stride = stride;
  cr->elt_size = elt_size;
  cr->comp_size = comp_size;
  cr->n_elts[0] = n_elts;
  cr->n_elts[1] = 0;

  cr->comm = comm;
  cr->rank_id = rank_id;
  cr->n_ranks = n_ranks;

  /* Compute data size and alignment */

  cr->dest_id_shift = 2*sizeof(int);

  if (cr->dest_id_shift % sizeof(cs_gnum_t) != 0)
    cr->dest_id_shift
      += (sizeof(cs_gnum_t) - (cr->dest_id_shift % sizeof(cs_gnum_t)));

  if (cr->dest_id_datatype == CS_GNUM_TYPE) {
    cr->src_id_shift = cr->dest_id_shift + sizeof(cs_gnum_t);
    align_size = sizeof(cs_gnum_t);
  }
  else if (cr->dest_id_datatype == CS_LNUM_TYPE)
    cr->src_id_shift = cr->dest_id_shift + sizeof(cs_lnum_t);
  else
    cr->src_id_shift = cr->dest_id_shift;

  if (cr->add_src_id)
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
 * Create a crystal router for strided data.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   dest_id          <-- element destination id, global number, or NULL
 *   dest_rank        <-- destination rank for each element
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _crystal_router_t *
_crystal_create_s(size_t           n_elts,
                  int              stride,
                  cs_datatype_t    datatype,
                  cs_datatype_t    dest_id_datatype,
                  bool             add_src_id,
                  void            *elt,
                  void            *dest_id,
                  const int        dest_rank[],
                  MPI_Comm         comm)
{
  size_t i, j;
  unsigned char *_elt = elt;
  unsigned char *_dest_id = dest_id;

  /* Allocate structure */

  _crystal_router_t *cr = _crystal_create_meta_s(n_elts,
                                                 stride,
                                                 datatype,
                                                 dest_id_datatype,
                                                 add_src_id,
                                                 comm);
  /* Copy data */

  for (i = 0; i < n_elts; i++) {

    int *pr = (int *)(cr->buffer[0] + i*cr->comp_size);
    unsigned char *pe = cr->buffer[0] + i*cr->comp_size + cr->elt_shift;
    unsigned char *_psrc = _elt + i*cr->elt_size;

    pr[0] = dest_rank[i];
    pr[1] = cr->rank_id;
    for (j = 0; j < cr->elt_size; j++)
      pe[j] = _psrc[j];

  }

  /* Add metadata */

  if (cr->dest_id_datatype != CS_DATATYPE_NULL) {
    const size_t id_size =   (cr->dest_id_datatype == CS_GNUM_TYPE)
                           ? sizeof(cs_gnum_t) : sizeof(cs_lnum_t);
    for (i = 0; i < n_elts; i++) {
      unsigned char *pi = cr->buffer[0] + cr->comp_size + cr->dest_id_shift;
      unsigned char *_p_dest_id = _dest_id + i*id_size;
      for (j = 0; j < id_size; j++)
        pi[j] = _p_dest_id[j];
    }
  }

  if (cr->add_src_id) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (i = 0; i < n_elts; i++) {
      cs_lnum_t src_id = i;
      unsigned char *_src_id = (unsigned char *)(&src_id);
      unsigned char *pi = cr->buffer[0] + i*cr->comp_size + cr->src_id_shift;
      for (j = 0; j < id_size; j++)
        pi[j] = _src_id[j];
    }
  }

  return cr;
}

/*----------------------------------------------------------------------------
 * Create a crystal router for strided data using block information.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   dest_id_datatype <-- type of destination id (CS_GNUM_TYPE, CS_LNUM_TYPE
 *                        or CS_DATATYPE_NULL)
 *   add_src_id       <-- add source id metadata (id in elt array)
 *   elt              <-- element values
 *   elt_gnum         <-- global element numbers
 *   bi               <-- destination block distribution info
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static _crystal_router_t *
_crystal_create_from_block_s(size_t                 n_elts,
                             int                    stride,
                             cs_datatype_t          datatype,
                             cs_datatype_t          dest_id_datatype,
                             bool                   add_src_id,
                             void                  *elt,
                             const cs_gnum_t       *elt_gnum,
                             cs_block_dist_info_t   bi,
                             MPI_Comm               comm)
{
  size_t i, j;
  unsigned char *_elt = elt;

  /* Allocate structure */

  _crystal_router_t *cr = _crystal_create_meta_s(n_elts,
                                                 stride,
                                                 datatype,
                                                 dest_id_datatype,
                                                 add_src_id,
                                                 comm);
  /* Copy data */

  const int rank_step = bi.rank_step;
  const cs_gnum_t block_size = bi.block_size;

  for (i = 0; i < n_elts; i++) {

    int *pr = (int *)(cr->buffer[0] + i*cr->comp_size);
    unsigned char *pe = cr->buffer[0] + i*cr->comp_size + cr->elt_shift;
    unsigned char *_psrc = _elt + i*cr->elt_size;

    cs_gnum_t g_elt_id = elt_gnum[i] -1;
    cs_gnum_t _dest_rank = g_elt_id / block_size;
    int dest_rank = _dest_rank*rank_step;

    pr[0] = dest_rank;
    pr[1] = cr->rank_id;
    for (j = 0; j < cr->elt_size; j++)
      pe[j] = _psrc[j];

  }

  /* Add metadata */

  if (cr->dest_id_datatype == CS_GNUM_TYPE) {
    const size_t gnum_size = sizeof(cs_gnum_t);
    for (i = 0; i < n_elts; i++) {
      cs_gnum_t g_elt_num = elt_gnum[j];
      unsigned char *_g_elt_num = (unsigned char *)(&g_elt_num);
      unsigned char *pi = cr->buffer[0] + cr->comp_size + cr->dest_id_shift;
      for (j = 0; j < gnum_size; j++)
        pi[j] = _g_elt_num[j];
    }
  }

  else if (cr->dest_id_datatype == CS_LNUM_TYPE) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (i = 0; i < n_elts; i++) {
      cs_gnum_t g_elt_id = elt_gnum[i] -1;
      cs_lnum_t b_elt_id = g_elt_id % block_size;;
      unsigned char *_b_elt_id = (unsigned char *)(&b_elt_id);
      unsigned char *pi = cr->buffer[0] + cr->comp_size + cr->dest_id_shift;
      for (j = 0; j < id_size; j++)
        pi[j] = _b_elt_id[j];
    }
  }

  if (cr->add_src_id) {
    const size_t id_size = sizeof(cs_lnum_t);
    for (i = 0; i < n_elts; i++) {
      cs_lnum_t src_id = j;
      unsigned char *_src_id = (unsigned char *)(&src_id);
      unsigned char *pi = cr->buffer[0] + i*cr->comp_size + cr->src_id_shift;
      for (j = 0; j < id_size; j++)
        pi[j] = _src_id[j];
    }
  }

  return cr;
}

/*----------------------------------------------------------------------------
 * Destroy a crystal router.
 *
 * parameters:
 *   cr <-> pointer to pointer to crystal router structure
 *---------------------------------------------------------------------------*/

static void
_crystal_destroy(_crystal_router_t **cr)
{
  if (cr != NULL) {
    _crystal_router_t *_cr = *cr;
    if (_cr->comp_type != MPI_BYTE)
      MPI_Type_free(&(_cr->comp_type));
    BFT_FREE(_cr->buffer[1]);
    BFT_FREE(_cr->buffer[0]);
    BFT_FREE(_cr);
  }
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
_crystal_partition_strided(_crystal_router_t  *cr,
                           int                 send_id,
                           int                 cutoff)
{
  cs_lnum_t i;

  cs_lnum_t n0 = 0, n1 = 0;

  const cs_lnum_t n = cr->n_elts[0];
  const int id0 = (send_id + 1) % 2;;
  const int id1 = send_id;
  const size_t comp_size = cr->comp_size;

  assert(send_id == 0 || send_id == 1);

  if (cr->buffer_size[1] < cr->buffer_size[0]) {
    cr->buffer_size[1] = cr->buffer_size[0];
    BFT_REALLOC(cr->buffer[1], cr->buffer_size[1], unsigned char);
  }

  for (i = 0; i < n; i++) {
    unsigned char *src = cr->buffer[0] + i*comp_size;
    int *r = (int *)src;
    if (r[0] < cutoff) {
      unsigned char *dest = cr->buffer[id0] + n0*comp_size;
      memmove(dest, src, comp_size);
      n0++;
    }
    else {
      unsigned char *dest = cr->buffer[id1] + n1*comp_size;
      memmove(dest, src, comp_size);
      n1++;
    }
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
_crystal_sendrecv_strided(_crystal_router_t  *cr,
                          int                 target,
                          int                 n_recv)
{
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
              "  Message to send would has size too large for C int: %llu",
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
_crystal_exchange_strided(_crystal_router_t  *cr)
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
_crystal_sort_by_source_rank_strided(_crystal_router_t  *cr)
{
  const cs_lnum_t n = cr->n_elts[0];

  if (n > 0)
    qsort(cr->buffer[0], n, cr->comp_size, &_compare_strided);
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor for strided data.
 *
 * \param[in]  n_elts     number of elements
 * \param[in]  stride     number of values per entity (interlaced)
 * \param[in]  datatype   type of data considered
 * \param[in]  elt        element values
 * \param[in]  dest_rank  destination rank for each element
 * \param[in]  comm       associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_s(size_t          n_elts,
                       int             stride,
                       cs_datatype_t   datatype,
                       void           *elt,
                       const int       dest_rank[],
                       MPI_Comm        comm)
{
  cs_all_to_all_t *d;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_all_to_all_calls[0] == 0) {
    int i;
    int n_timers = sizeof(_all_to_all_timers)/sizeof(_all_to_all_timers[0]);
    for (i = 0; i < n_timers; i++)
      CS_TIMER_COUNTER_INIT(_all_to_all_timers[i]);
  }

  /* Allocate structure */

  BFT_MALLOC(d, 1, cs_all_to_all_t);

  /* Create associated sub-structure */

  d->strided = true;
  d->cr = NULL;
  d->dc = NULL;

  if (_all_to_all_type == CS_ALL_TO_ALL_CRYSTAL_ROUTER)
    d->cr = _crystal_create_s(n_elts,
                              stride,
                              datatype,
                              CS_DATATYPE_NULL,
                              false,
                              elt,
                              NULL,
                              dest_rank,
                              comm);

  else
    d->dc = _alltoall_caller_create_s(n_elts,
                                      stride,
                                      datatype,
                                      CS_DATATYPE_NULL,
                                      false,
                                      elt,
                                      NULL,
                                      dest_rank,
                                      comm);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers, &t0, &t1);
  _all_to_all_calls[0] += 1;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor for strided data with additional
 *  metadata.
 *
 * This variant allows optional tracking of destination ids or global
 * numbers associated with elements, as well as their source ids.
 *
 * In cases where those arrays are required and already available, this
 * may avoid the need for a specific element values buffer mixing actual
 * data values and numbering metadata. It also makes extraction of the
 * metadata easier using \ref cs_all_to_all_get_id_pointers.
 *
 * \param[in]  n_elts            number of elements
 * \param[in]  stride            number of values per entity (interlaced)
 * \param[in]  datatype          type of data considered
 * \param[in]  dest_id_datatype  type of destination id (CS_GNUM_TYPE,
 *                               CS_LNUM_TYPE or CS_DATATYPE_NULL depending
 *                               on elt_id values)
 * \param[in]  add_src_id        add source id metadata (id in elt array)
 * \param[in]  elt               element values
 * \param[in]  dest_id           element destination id, global number, or NULL
 * \param[in]  dest_rank         destination rank for each element
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_with_ids_s(size_t            n_elts,
                                int               stride,
                                cs_datatype_t     datatype,
                                cs_datatype_t     dest_id_datatype,
                                bool              add_src_id,
                                void             *elt,
                                void             *dest_id,
                                const int         dest_rank[],
                                MPI_Comm          comm)
{
  cs_all_to_all_t *d;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_all_to_all_calls[0] == 0) {
    int i;
    int n_timers = sizeof(_all_to_all_timers)/sizeof(_all_to_all_timers[0]);
    for (i = 0; i < n_timers; i++)
      CS_TIMER_COUNTER_INIT(_all_to_all_timers[i]);
  }

  /* Allocate structure */

  BFT_MALLOC(d, 1, cs_all_to_all_t);

  /* Create associated sub-structure */

  d->strided = true;
  d->cr = NULL;
  d->dc = NULL;

  if (_all_to_all_type == CS_ALL_TO_ALL_CRYSTAL_ROUTER)
    d->cr = _crystal_create_s(n_elts,
                              stride,
                              datatype,
                              dest_id_datatype,
                              add_src_id,
                              elt,
                              dest_id,
                              dest_rank,
                              comm);

  else
    d->dc = _alltoall_caller_create_s(n_elts,
                                      stride,
                                      datatype,
                                      dest_id_datatype,
                                      add_src_id,
                                      elt,
                                      dest_id,
                                      dest_rank,
                                      comm);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers, &t0, &t1);
  _all_to_all_calls[0] += 1;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an all-to-all distributor for strided data with additional
 * metadata, with destination rank determined from global numbers and block
 * distribution information.
 *
 * This variant allows optional tracking of destination ids or global
 * numbers associated with elements, as well as their source ids.
 *
 * In cases where those arrays are required and already available, this
 * may avoid the need for a specific element values buffer mixing actual
 * data values and numbering metadata. It also makes extraction of the
 * metadata easier using \ref cs_all_to_all_get_id_pointers.
 *
 * \param[in]  n_elts            number of elements
 * \param[in]  stride            number of values per entity (interlaced)
 * \param[in]  datatype          type of data considered
 * \param[in]  dest_id_datatype  type of destination id (CS_GNUM_TYPE,
 *                               CS_LNUM_TYPE or CS_DATATYPE_NULL depending
 *                               on elt_id values)
 * \param[in]  add_src_id        add source id metadata (id in elt array)
 * \param[in]  elt               element values
 * \param[in]  elt_gnum          global element numbers
 * \param[in]  bi                destination block distribution info
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_all_to_all_create_from_block_s(size_t                 n_elts,
                                  int                    stride,
                                  cs_datatype_t          datatype,
                                  cs_datatype_t          dest_id_datatype,
                                  bool                   add_src_id,
                                  void                  *elt,
                                  const cs_gnum_t       *elt_gnum,
                                  cs_block_dist_info_t   bi,
                                  MPI_Comm               comm)
{
  cs_all_to_all_t *d;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Initialize timers if required */

  if (_all_to_all_calls[0] == 0) {
    int i;
    int n_timers = sizeof(_all_to_all_timers)/sizeof(_all_to_all_timers[0]);
    for (i = 0; i < n_timers; i++)
      CS_TIMER_COUNTER_INIT(_all_to_all_timers[i]);
  }

  /* Allocate structure */

  BFT_MALLOC(d, 1, cs_all_to_all_t);

  /* Create associated sub-structure */

  d->strided = true;

  if (_all_to_all_type == CS_ALL_TO_ALL_CRYSTAL_ROUTER)
    d->cr = _crystal_create_from_block_s(n_elts,
                                         stride,
                                         datatype,
                                         dest_id_datatype,
                                         add_src_id,
                                         elt,
                                         elt_gnum,
                                         bi,
                                         comm);

  else
    d->dc = _alltoall_caller_create_from_block_s(n_elts,
                                                 stride,
                                                 datatype,
                                                 dest_id_datatype,
                                                 add_src_id,
                                                 elt,
                                                 elt_gnum,
                                                 bi,
                                                 comm);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers, &t0, &t1);
  _all_to_all_calls[0] += 1;

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
      _crystal_destroy(&(_d->cr));
    else if (_d->dc != NULL) {
      _alltoall_caller_destroy(&(_d->dc));
    }

    BFT_FREE(_d);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(_all_to_all_timers, &t0, &t1);

    /* no increment to _all_to_all_calls[0] as create/destroy are grouped */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange data with an all-to-all distributor.
 *
 * Order of data from a same source rank is preserved.
 *
 * \param[in, out]  d   pointer to associated all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_exchange(cs_all_to_all_t  *d)
{
  if (d == NULL)
    return;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  if (d->cr != NULL) {
    if (d->strided)
      _crystal_exchange_strided(d->cr);
  }
  else if (d->dc != NULL) {
    if (d->strided)
      _alltoall_caller_exchange_strided(d->dc);
  }

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers+1, &t0, &t1);
  _all_to_all_calls[1] += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sort stride crystal router data by source rank.
 *
 * \param[in, out]  d   pointer to associated all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_sort_by_source_rank(cs_all_to_all_t  *d)
{
  if (d->cr != NULL) {

    cs_timer_t t0, t1;

    t0 = cs_timer_time();

    if (d->strided)
      _crystal_sort_by_source_rank_strided(d->cr);

    t1 = cs_timer_time();
    cs_timer_counter_add_diff(_all_to_all_timers+1, &t0, &t1);
    _all_to_all_calls[3] += 1;
  }
  /* For MPI_Alltoall/MPI_Alltoallv, data is already sorted */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements associated with all-to-all distributor.
 *
 * The number of elements is the number of elements received after exchange.
 *
 * \param[in]  d   pointer to associated all-to-all distributor
 *
 * \return  number of elements associated with distributor.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_all_to_all_n_elts(const cs_all_to_all_t  *d)
{
  cs_lnum_t retval = 0;

  if (d != NULL) {
    if (d->cr != NULL)
      retval = d->cr->n_elts[0];
    else if (d->dc != NULL)
      retval = d->dc->recv_size;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Swap source and destination ranks of all-to-all distributor.
 *
 * \param[in, out]  d   pointer to associated all-to-all distributor
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_swap_src_dest(cs_all_to_all_t  *d)
{
  if (d == NULL)
    return;

  cs_timer_t t0, t1;

  t0 = cs_timer_time();

  /* Case for Crystal router */

  if (d->cr != NULL) {

    size_t i;
    size_t comp_size = d->cr->comp_size;
    unsigned char *p = d->cr->buffer[0];

    for (i = 0; i < d->cr->n_elts[0]; i++) {
      unsigned char *p1 = p + i*comp_size;
      int *r = (int *)p1;
      int t = r[0];
      r[0] = r[1];
      r[1] = t;
    }

  }

  /* Case for MPI_Alltoall/MPI_Alltoallv */

  else if (d->dc != NULL) {

    size_t tmp_size[2] = {d->dc->send_size, d->dc->recv_size};
    unsigned char  *tmp_buffer = d->dc->recv_buffer;
    int *tmp_count = d->dc->recv_count;
    int *tmp_displ = d->dc->recv_displ;
    int *tmp_rank = d->dc->src_rank;

    d->dc->send_size = tmp_size[1];
    d->dc->recv_size = tmp_size[0];

    d->dc->recv_buffer = d->dc->send_buffer;
    d->dc->recv_count = d->dc->send_count;
    d->dc->recv_displ = d->dc->send_displ;
    d->dc->src_rank = d->dc->dest_rank;

    d->dc->send_buffer = tmp_buffer;
    d->dc->send_count = tmp_count;
    d->dc->send_displ = tmp_displ;
    d->dc->dest_rank = tmp_rank;

  }

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(_all_to_all_timers+1, &t0, &t1);
  _all_to_all_calls[2] += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to data elements associated with an all-to-all
 * distributor.
 *
 * This allows modification and/or extraction of those elements.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * \param[in]   d            pointer to associated all-to-all distributor
 * \param[out]  data_stride  stride (in bytes) between data items
 * \param[out]  data         pointer to data items
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_data_pointer(cs_all_to_all_t   *d,
                               size_t            *data_stride,
                               unsigned char    **data)
{
  *data_stride = 0;
  *data  = NULL;

  if (d == NULL)
    return;

  if (d->strided) {

    /* Case for Crystal router */

    if (d->cr != NULL) {
      *data_stride = d->cr->comp_size;
      *data = d->cr->buffer[0] + d->cr->elt_shift;
    }

    /* Case for MPI_Alltoall/MPI_Alltoallv */

    else if (d->dc != NULL) {
      *data_stride = d->dc->comp_size;
      *data = d->dc->recv_buffer + d->dc->elt_shift;
    }

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              "%s is only available for strided (not indexed) data.",
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to ranks of elements associated with an all-to-all
 *  distributor.
 *
 * This allows modification and/or extraction of those ranks.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * \param[in]   d            pointer to associated all-to-all distributor
 * \param[out]  rank_stride  stride (in integers) between rank values
 * \param[out]  src_rank     pointer to source rank values (or NULL)
 * \param[out]  dest_rank    pointer to destination rank values (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_rank_pointers(cs_all_to_all_t   *d,
                                size_t            *rank_stride,
                                int              **src_rank,
                                int              **dest_rank)
{
  *rank_stride = 0;
  if (src_rank != NULL)
    *src_rank = NULL;
  if (dest_rank != NULL)
    *dest_rank = NULL;

  if (d == NULL)
    return;

  if (d->strided) {

    /* Case for Crystal router */

    if (d->cr != NULL) {
      const size_t comp_size = d->cr->comp_size;
      unsigned char *buffer = d->cr->buffer[0];
      if (rank_stride != NULL)
        *rank_stride = comp_size / sizeof(int);
      if (src_rank != NULL)
        *src_rank = (int *)(buffer + sizeof(int));
      if (dest_rank != NULL)
        *dest_rank = (int *)buffer;
    }

    /* Case for MPI_Alltoall/MPI_Alltoallv (build arrays upon request) */

    else if (d->dc != NULL) {

      int i;
      cs_lnum_t j;
      _mpi_all_to_all_caller_t *dc = d->dc;

      if (rank_stride != NULL)
        *rank_stride = 1;

      if (src_rank != NULL) {
        if (dc->src_rank == NULL) {
          BFT_MALLOC(dc->src_rank, dc->recv_size, int);
          for (i = 0; i < dc->n_ranks; i++) {
            for (j = dc->recv_displ[i]; j < dc->recv_displ[i+1]; j++)
              dc->src_rank[j] = i;
          }
        }
        *src_rank = dc->src_rank;
      }

      if (dest_rank != NULL) {
        if (dc->dest_rank == NULL) {
          BFT_MALLOC(dc->dest_rank, dc->send_size, int);
          for (i = 0; i < dc->n_ranks; i++) {
            for (j = dc->send_displ[i]; j < dc->send_displ[i+1]; j++)
              dc->dest_rank[j] = i;
          }
        }
        *dest_rank = dc->dest_rank;
      }

    }

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              "%s is only available for strided (not indexed) data.",
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to source or destination rank element ids associated
 *  with an all-to-all distributor.
 *
 * If a requested type of id is not available (depending on the all-to-all
 * distributor creation function and options), the matching pointer will
 * be set to NULL.
 *
 * This allows modification and/or extraction of those ids, though it is
 * intended primarily for identification.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * \param[in]   d          pointer to associated all-to-all distributor
 * \param[out]  id_stride  stride (in integers) between id items
 * \param[out]  dest_id    pointer to destination ids (or NULL)
 * \param[out]  src_id     pointer to source ids (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_id_pointers(cs_all_to_all_t   *d,
                              size_t            *id_stride,
                              cs_lnum_t        **dest_id,
                              cs_lnum_t        **src_id)
{
  *id_stride = 0;
  if (dest_id != NULL)
    *dest_id = NULL;
  if (src_id != NULL)
    *src_id = NULL;

  if (d == NULL)
    return;

  if (d->strided) {

    /* Case for Crystal router */

    if (d->cr != NULL) {
      *id_stride = d->cr->comp_size / sizeof(cs_lnum_t);
      if (dest_id != NULL) {
        if (d->cr->dest_id_datatype == CS_LNUM_TYPE)
          *dest_id = (cs_lnum_t *)(d->cr->buffer[0] + d->cr->dest_id_shift);
        else
          *dest_id = NULL;
      }
      if (src_id != NULL) {
        if (d->cr->add_src_id)
          *src_id = (cs_lnum_t *)(d->cr->buffer[0] + d->cr->src_id_shift);
        else
          *src_id = NULL;
      }
    }

    /* Case for MPI_Alltoall/MPI_Alltoallv */

    else if (d->dc != NULL) {
      *id_stride = d->dc->comp_size / sizeof(cs_lnum_t);
      if (dest_id != NULL) {
        if (d->dc->dest_id_datatype == CS_LNUM_TYPE)
          *dest_id = (cs_lnum_t *)(d->dc->recv_buffer);
        else
          *dest_id = NULL;
      }
      if (src_id != NULL) {
        if (d->dc->add_src_id)
          *src_id = (cs_lnum_t *)(d->dc->recv_buffer + d->dc->src_id_shift);
        else
          *src_id = NULL;
      }
    }
  }

  else
    bft_error(__FILE__, __LINE__, 0,
              "%s is only available for strided (not indexed) data.",
              __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to element global numbers associated  with an all-to-all
 *  distributor.
 *
 * If this data is not available (depending on the all-to-all distributor
 * creation function and options), the matching pointer will
 * be set to NULL.
 *
 * This allows modification and/or extraction of those numbers.
 *
 * Note that depending on the distributor type used, the rank metadata
 * and data may be interleaved, so the corresponding pointers point
 * to strided, interleaved data.
 *
 * \param[in]   d            pointer to associated all-to-all distributor
 * \param[out]  gnum_stride  stride (in integers) between element global numbers
 * \param[out]  gnum         pointer to global numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_all_to_all_get_gnum_pointer(cs_all_to_all_t   *d,
                               size_t            *gnum_stride,
                               cs_gnum_t        **gnum)
{
  *gnum_stride = 0;
  *gnum = NULL;

  if (d == NULL)
    return;

  if (d->strided) {

    /* Case for Crystal router */

    if (d->cr != NULL) {
      *gnum_stride = d->cr->comp_size / sizeof(cs_gnum_t);
      if (d->cr->dest_id_datatype == CS_GNUM_TYPE)
        *gnum = (cs_gnum_t *)(d->cr->buffer[0] + d->cr->dest_id_shift);
      else
        *gnum = NULL;
    }

    /* Case for MPI_Alltoall/MPI_Alltoallv */

    else if (d->dc != NULL) {
      *gnum_stride = d->dc->comp_size / sizeof(cs_gnum_t);
      if (d->dc->dest_id_datatype == CS_GNUM_TYPE)
        *gnum = (cs_gnum_t *)(d->dc->recv_buffer);
    }

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              "%s is only available for strided (not indexed) data.",
              __func__);
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
#if defined(HAVE_MPI)

  int i;
  size_t name_width = 0;

  const char *method_name[] = {N_("MPI_Alltoall and MPI_Alltoallv"),
                               N_("Crystal Router algorithm")};
  const char *operation_name[] = {N_("Construction/destruction:"),
                                  N_("Exchange:"),
                                  N_("Swap source and destination:"),
                                  N_("Sort by source rank:"),
                                  N_("Copy exchanged data:")};

  if (_all_to_all_calls[0] <= 0)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nInstrumented all-to-all operations (using %s):\n\n"),
                _(method_name[_all_to_all_type]));

  /* Determine width */

  for (i = 0; i < 5; i++) {
    if (_all_to_all_calls[i] > 0) {
      size_t l = cs_log_strlen(_(operation_name[i]));
      name_width = CS_MAX(name_width, l);
    }
  }
  name_width = CS_MIN(name_width, 63);

  /* Print times */

  for (i = 0; i < 5; i++) {

    if (_all_to_all_calls[i] > 0) {

      char tmp_s[64];
      double wtimes = (_all_to_all_timers[i]).wall_nsec*1e-9;

      cs_log_strpad(tmp_s, _(operation_name[i]), name_width, 64);

      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("  %s %12.5f s, %lu calls\n"),
                    tmp_s, wtimes, (unsigned long)(_all_to_all_calls[i]));

    }

  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

#endif /* defined(HAVE_MPI) */
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
