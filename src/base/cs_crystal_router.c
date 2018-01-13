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

#define _CR_DEBUG_DUMP 0

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
  size_t          n_vals_shift;      /* starting byte for element count
                                        (for indexed cases) */
  size_t          elt_shift;         /* starting byte for element data */

  size_t          elt_size;          /* element size */
  size_t          comp_size;         /* composite metadata + element size if
                                        strided, metadata size otherwise */
  size_t          n_elts[2];         /* number of elements in partition */
  size_t          n_vals[2];         /* number of data values in partition */
  size_t          buffer_size[2];    /* buffer size values */
  unsigned char  *buffer[2];

  MPI_Comm        comm;              /* associated MPI communicator */
  MPI_Datatype    mpi_type;          /* Associated MPI datatype */
  size_t          mpi_type_size;     /* Associated MPI datatype size */
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
 * Return Crystal Router datatype count matching a given number of values.
 *
 * parameters:
 *   cr      <-- associated crystal router structure
 *   n_elts  <-- number of elements
 *   n_vals  <-- number of associated values (used only for indexed cases)
 *
 * returns:
 *   count of associated datatype elements for matching buffer size
 *---------------------------------------------------------------------------*/

static  size_t
_comm_type_count(cs_crystal_router_t  *cr,
                 cs_lnum_t             n_elts,
                 cs_lnum_t             n_vals)
{
  size_t retval;

  if (cr->n_vals_shift == 0) { /* strided */
    retval = n_elts;
  }
  else { /* indexed */
    size_t n = n_elts*cr->comp_size + n_vals*cr->elt_size;
    retval = n / cr->mpi_type_size;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return number of values matching a given number of elements and MPI
 * datatype values.
 *
 * parameters:
 *   cr         <-- associated crystal router structure
 *   n_elts     <-- number of elements
 *   n_mpi_vals <-- number of MPI element
 *
 * returns:
 *   count of associated datatype elements for matching buffer size
 *---------------------------------------------------------------------------*/

static  size_t
_comm_type_count_to_n_vals(cs_crystal_router_t  *cr,
                           cs_lnum_t             n_elts,
                           int                   n_mpi_elts)
{
  size_t retval;

  if (cr->n_vals_shift == 0) { /* strided */
    retval = n_elts * cr->stride;
  }
  else { /* indexed */
    size_t n = n_mpi_elts*cr->mpi_type_size - n_elts*cr->comp_size;
    retval = n / cr->elt_size;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return Crystal Router datatype count matching a given number of values.
 *
 * parameters:
 *   cr      <-- associated crystal router structure
 *   n_elts  <-- number of elements
 *   n_vals  <-- number of associated values (used only for indexed cases)
 *
 * returns:
 *   count of associated datatype elements for matching buffer size
 *---------------------------------------------------------------------------*/

static inline size_t
_data_size(cs_crystal_router_t  *cr,
           cs_lnum_t             n_elts,
           cs_lnum_t             n_vals)
{
  size_t retval;

  if (cr->n_vals_shift == 0) /* strided */
    retval = n_elts * cr->comp_size;
  else /* indexed */
    retval = n_elts*cr->comp_size + n_vals*cr->elt_size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Allocation and base settings for a crystal router.
 *
 * parameters:
 *   n_elts           <-- number of elements
 *   stride           <-- number of values per entity (interlaced)
 *   datatype         <-- type of data considered
 *   flags            <-- ordering and metadata flags
 *   comm             <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static cs_crystal_router_t *
_crystal_create(size_t         n_elts,
                int            flags,
                MPI_Comm       comm)
{
  int rank_id, n_ranks;
  cs_crystal_router_t *cr = NULL;

  const size_t align_size = sizeof(cs_lnum_t);

  /* Communicator info */

  MPI_Comm_rank(comm, &rank_id);
  MPI_Comm_size(comm, &n_ranks);

  /* Allocate structure */

  BFT_MALLOC(cr, 1, cs_crystal_router_t);

  cr->flags = flags;

  cr->stride = 0;
  cr->elt_size = 0;
  cr->comp_size = 0;
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

  cr->n_vals_shift = 0;

  for (int i = 0; i < 2; i++) {
    cr->n_vals[i] = 0;
    cr->buffer_size[i] = 0;
    cr->buffer[i] = NULL;
  }

  return cr;
}

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
  /* Allocate structure */

  cs_crystal_router_t *cr = _crystal_create(n_elts, flags, comm);

  size_t elt_size = cs_datatype_size[datatype]*stride;
  size_t align_size = sizeof(cs_lnum_t);

  cr->datatype = (stride > 0) ? datatype : CS_DATATYPE_NULL;

  cr->stride = (stride > 0) ? stride : 1;
  cr->elt_size = elt_size;

  /* Compute data size and alignment */

  cr->comp_size = cr->elt_shift + elt_size;

  if (elt_size % align_size)
    cr->comp_size += align_size - (elt_size % align_size);;

  /* Create associated MPI datatype */

  cr->mpi_type_size = cr->comp_size;

  MPI_Type_contiguous(cr->mpi_type_size, MPI_BYTE, &(cr->mpi_type));
  MPI_Type_commit(&(cr->mpi_type));

  /* Allocate buffers */

  cr->buffer_size[0] = n_elts*cr->comp_size;
  cr->buffer_size[1] = 0;
  cr->buffer_size[1] = 0;
  BFT_MALLOC(cr->buffer[0], cr->buffer_size[0], unsigned char);
  memset(cr->buffer[0], 0, cr->buffer_size[0]);
  cr->buffer[1] = NULL;

  return cr;
}

/*----------------------------------------------------------------------------
 * First stage of creation for a crystal router for indexed data.
 *
 * parameters:
 *   n_elts     <-- number of elements
 *   datatype   <-- type of data considered
 *   elt_idx    <-- element values start and past-the-last index
 *   flags      <-- ordering and metadata flags
 *   comm       <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

static cs_crystal_router_t *
_crystal_create_meta_i(size_t            n_elts,
                       cs_datatype_t     datatype,
                       const cs_lnum_t  *elt_idx,
                       int               flags,
                       MPI_Comm          comm)
{
  assert(datatype != CS_DATATYPE_NULL);

  /* Allocate structure */

  cs_crystal_router_t *cr = _crystal_create(n_elts, flags, comm);

  size_t comp_size = 0;
  size_t elt_size = cs_datatype_size[datatype];
  size_t align_size = sizeof(cs_lnum_t);

  /* Ensure alignement on integer size at least */

  if (elt_size % sizeof(int) > 0)
    comp_size += sizeof(int) - elt_size % sizeof(int);

  cr->datatype = datatype;

  cr->elt_size = elt_size;

  /* Compute data size and alignment */

  cr->n_vals_shift = cr->elt_shift;
  cr->elt_shift = cr->n_vals_shift + cs_datatype_size[CS_LNUM_TYPE];

  if (cr->elt_shift < cr->elt_size)
    cr->elt_shift = cr->elt_size;

  if (cr->elt_shift % align_size)
    cr->elt_shift += align_size - (cr->elt_shift % align_size);

  cr->comp_size = cr->elt_shift;

  /* Create associated MPI datatype */

  cr->mpi_type_size = CS_MIN(cr->comp_size, cr->elt_size);
  while (cr->comp_size % cr->mpi_type_size || cr->elt_size % cr->mpi_type_size)
    cr->mpi_type_size--;

  MPI_Type_contiguous(cr->mpi_type_size, MPI_BYTE, &(cr->mpi_type));
  MPI_Type_commit(&(cr->mpi_type));

  /* Allocate buffers */

  cr->buffer_size[0] = n_elts*cr->comp_size + elt_idx[n_elts]*elt_size;
  cr->buffer_size[1] = 0;
  cr->n_vals[0] = elt_idx[n_elts];
  cr->n_vals[1] = 0;
  BFT_MALLOC(cr->buffer[0], cr->buffer_size[0], unsigned char);
  memset(cr->buffer[0], 0, cr->buffer_size[0]);
  cr->buffer[1] = NULL;

  return cr;
}

#if _CR_DEBUG_DUMP /* _dump functions for debugging */

/*----------------------------------------------------------------------------
 * Dump data element associated with a Crystal Router.
 *
 * parameters:
 *   <-- cr pointer to associated Crystal Router
 *   <-- e  pointer to element
 *----------------------------------------------------------------------------*/

static void
_dump_e(const cs_crystal_router_t  *cr,
        const void                 *e)
{
  switch (cr->datatype) {

  case CS_FLOAT:
    {
      const float *v = (const float *)e;
      bft_printf(" %g", (double)(*v));
    }
    break;
  case CS_DOUBLE:
    {
      const double *v = (const double *)e;
      bft_printf(" %g", (*v));
    }
    break;
  case CS_LNUM_TYPE:
    {
      const cs_lnum_t *v = (const cs_lnum_t *)e;
      bft_printf(" %d", (int)(*v));
    }
    break;
  case CS_GNUM_TYPE:
    {
      const cs_gnum_t *v = (const cs_gnum_t *)e;
      bft_printf(" %llu", (unsigned long long)(*v));
    }
    break;
  default:
    {
      const unsigned char *pe = (const unsigned char *)e;
      size_t elt_size = cs_datatype_size[cr->datatype];
      if (cr->elt_size > 0)
        bft_printf(" %x", (int)pe[0]);
      for (size_t k = 1; k < elt_size; k++)
        bft_printf(":%x", (int)pe[k]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Dump data elements associated with a strided Crystal Router.
 *
 * parameters:
 *   <-- cr       pointer to associated Crystal Router
 *   <-- comment  associated comment
 *----------------------------------------------------------------------------*/

static void
_dump_s(const cs_crystal_router_t  *cr,
        const char                 *comment)
{
  if (comment != NULL)
    bft_printf("\n%s\n", comment);

  size_t elt_size = cr->elt_size;
  if (cr->stride > 0)
    elt_size = cr->elt_size/cr->stride;

  for (int ip = 0; ip < 2; ip++) {

    const size_t n_elts = cr->n_elts[ip];

    if (n_elts > 0)
      bft_printf("crystal router partition %d: %d elements "
                 "(stride %d, size %d)\n",
                 ip, (int)n_elts, (int)cr->stride, (int)elt_size);

    /* Extract data */

    for (size_t i = 0; i < n_elts; i++) {

      unsigned const char *p_s = cr->buffer[ip] + i*cr->comp_size;

      const int *cr_dest_rank = (const int *)(p_s);

      bft_printf("  %d\n"
                 "    dest_rank:  %d\n", (int)i, *cr_dest_rank);

      if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_RANK) {
        const int *cr_src_rank = (const int *)(p_s + sizeof(int));
        bft_printf("    src_rank:   %d\n", (int)(*cr_src_rank));
      }

      if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
        const cs_lnum_t *cr_dest_id
          = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
        bft_printf("    dest_id:    %d\n", (int)(*cr_dest_id));
      }

      if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_ID) {
        const cs_lnum_t *cr_src_id
          = (const cs_lnum_t *)(p_s + cr->src_id_shift);
        bft_printf("    src_id:     %d\n", (int)(*cr_src_id));
      }

      /* Data */

      if (cr->elt_size > 0) {
        unsigned const char *pe = p_s + cr->elt_shift;
        bft_printf("    data      :");
        for (size_t j = 0; j < cr->stride; j++)
          _dump_e(cr, pe + elt_size*j);
      }
      bft_printf("\n");
      bft_printf_flush();

    }
  }
}

/*----------------------------------------------------------------------------
 * Dump data elements associated with an indexed Crystal Router.
 *
 * parameters:
 *   <-- cr       pointer to associated Crystal Router
 *   <-- comment  associated comment
 *----------------------------------------------------------------------------*/

static void
_dump_i(const cs_crystal_router_t  *cr,
        const char                 *comment)
{
  if (comment != NULL)
    bft_printf("\n%s\n", comment);

  for (int ip = 0; ip < 2; ip++) {

    const size_t n_elts = cr->n_elts[ip];

    if (n_elts > 0)
      bft_printf("crystal router partition %d: %d elements\n",
                 ip, (int)n_elts);

    unsigned const char *p_s = cr->buffer[ip];

    /* Extract data */

    cs_lnum_t  s_idx = 0;

    for (size_t i = 0; i < n_elts; i++) {

      const int *cr_dest_rank = (const int *)(p_s);

      bft_printf("  %d\n"
                 "    dest_rank:  %d\n", (int)i, *cr_dest_rank);

      if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_RANK) {
        const int *cr_src_rank = (const int *)(p_s + sizeof(int));
        bft_printf("    src_rank:   %d\n", (int)(*cr_src_rank));
      }

      if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
        const cs_lnum_t *cr_dest_id
          = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
        bft_printf("    dest_id:    %d\n", (int)(*cr_dest_id));
      }

      if (cr->flags & CS_CRYSTAL_ROUTER_ADD_SRC_ID) {
        const cs_lnum_t *cr_src_id
          = (const cs_lnum_t *)(p_s + cr->src_id_shift);
        bft_printf("    src_id:     %d\n", (int)(*cr_src_id));
      }

      /* Number of elements and optional data index */

      const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
      const cs_lnum_t n_sub = *pn;
      const size_t _sub_size = cr->elt_size*n_sub;

      bft_printf("    data index: %d %d\n"
                 "    data      : ", (int)s_idx, (int)(s_idx+n_sub));

      /* Data */

      unsigned const char *pe = p_s + cr->elt_shift;
      for (cs_lnum_t j = 0; j < n_sub; j++)
        _dump_e(cr, pe + cr->elt_size*j);
      bft_printf("\n");
      bft_printf_flush();

      s_idx += n_sub;
      p_s += cr->comp_size + _sub_size;

    }
  }
}

/*----------------------------------------------------------------------------
 * Dump data elements associated with a Crystal Router.
 *
 * parameters:
 *   <-- cr       pointer to associated Crystal Router
 *   <-- comment  associated comment
 *----------------------------------------------------------------------------*/

static void
_dump(const cs_crystal_router_t  *cr,
      const char                 *comment)
{
  if (cr->n_vals_shift == 0)
    _dump_s(cr, comment);
  else
    _dump_i(cr, comment);
}

#endif /* _dump functions for debugging */

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

  cr->n_vals[id0] = n0*cr->stride;
  cr->n_vals[id1] = n1*cr->stride;
}

/*----------------------------------------------------------------------------
 * Partition indexed data for exchange with a crystal router.
 *
 * parameters:
 *   cr      <-> associated crystal router structure
 *   send_id <-> id of "send" buffer" (1 if low = buf0, 0 if low = buf1)
 *   cutoff  <-- cutoff rank
 *---------------------------------------------------------------------------*/

static void
_crystal_partition_indexed(cs_crystal_router_t  *cr,
                           int                   send_id,
                           int                   cutoff)
{
  cs_lnum_t i;

  cs_lnum_t n0 = 0, n1 = 0;

  const cs_lnum_t n = cr->n_elts[0];
  const int id0 = (send_id + 1) % 2;
  const int id1 = send_id;
  const size_t comp_size = cr->comp_size;
  const size_t elt_size = cr->elt_size;
  const size_t n_vals_shift = cr->n_vals_shift;

  assert(send_id == 0 || send_id == 1);

  if (cr->buffer_size[1] < cr->buffer_size[0]) {
    cr->buffer_size[1] = cr->buffer_size[0];
    BFT_REALLOC(cr->buffer[1], cr->buffer_size[1], unsigned char);
  }

  unsigned char *src = cr->buffer[0];

  size_t r0_shift = 0;
  size_t r1_shift = 0;

  if (id0 == 0) {
    for (i = 0; i < n; i++) {
      int *r = (int *)src;
      cs_lnum_t *n_sub = (cs_lnum_t *)(src + n_vals_shift);
      size_t sub_size = comp_size + n_sub[0]*elt_size;
      if (r[0] >= cutoff)
        break;
      r0_shift += sub_size;
      src += sub_size;
    }
    n0 = i;
  }
  else {
    for (i = 0; i < n; i++) {
      int *r = (int *)src;
      cs_lnum_t *n_sub = (cs_lnum_t *)(src + n_vals_shift);
      size_t sub_size = comp_size + n_sub[0]*elt_size;
      if (r[0] < cutoff)
        break;
      r1_shift += sub_size;
      src += sub_size;
    }
    n1 = i;
  }

  while (i < n) {
    int *r = (int *)src;
    cs_lnum_t *n_sub = (cs_lnum_t *)(src + n_vals_shift);
    size_t sub_size = comp_size + n_sub[0]*elt_size;
    if (r[0] < cutoff) {
      unsigned char *dest = cr->buffer[id0] + r0_shift;
      for (size_t j = 0; j < sub_size; j++)
        dest[j] = src[j];
      r0_shift += sub_size;
      n0++;
    }
    else {
      unsigned char *dest = cr->buffer[id1] + r1_shift;
      for (size_t j = 0; j < sub_size; j++)
        dest[j] = src[j];
      r1_shift += sub_size;
      n1++;
    }
    i++;
    src += sub_size;
  }

  assert(n0 + n1 == n);

  cr->n_elts[id0] = n0;
  cr->n_elts[id1] = n1;

  cr->n_vals[id0] = (r0_shift - n0*comp_size)/elt_size;
  cr->n_vals[id1] = (r1_shift - n1*comp_size)/elt_size;
}

/*----------------------------------------------------------------------------
 * Send and receive data with a crystal router for one stage.
 *
 * parameters:
 *   cr        <-> associated crystal router structure
 *   target    <-- base target rank to send/receive from
 *   n_recv    <-- number of ranks to receive from
 *---------------------------------------------------------------------------*/

static  void
_crystal_sendrecv(cs_crystal_router_t  *cr,
                  int                   target,
                  int                   n_recv)
{
  cs_timer_t t0 = cs_timer_time();

  assert(n_recv <= 2);

  cs_lnum_t send_size[2];
  uint64_t test_size;

  MPI_Status status[3];
  MPI_Request request[3] = {MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL};
  cs_lnum_t recv_size[4] = {0, 0, 0, 0};

  /* Send message to target process */

  test_size = _comm_type_count(cr, cr->n_elts[1], cr->n_vals[1]);

  send_size[0] = cr->n_elts[1];
  send_size[1] = test_size;

  if ((uint64_t)send_size[1] != test_size)
    bft_error(__FILE__, __LINE__, 0,
              "Crystal router:"
              "  Message to send would have size too large for C int: %llu",
              (unsigned long long)test_size);

  MPI_Isend(&send_size, 2, CS_MPI_LNUM, target, cr->rank_id,
            cr->comm, &request[0]);

  for (int i = 0; i < n_recv; i++)
    MPI_Irecv(recv_size+i*2, 2, CS_MPI_LNUM, target+i, target+i,
              cr->comm, request+i+1);

  MPI_Waitall(n_recv + 1, request, status);

  size_t loc_size = _data_size(cr, cr->n_elts[0], cr->n_vals[0]);
  for (int i = 0; i < n_recv; i++)
    loc_size += _data_size(cr, recv_size[i*2], recv_size[i*2+1]);
  if (loc_size > cr->buffer_size[0]) {
    cr->buffer_size[0] = loc_size;
    BFT_REALLOC(cr->buffer[0], cr->buffer_size[0], unsigned char);
  }

  MPI_Isend(cr->buffer[1], send_size[1], cr->mpi_type,
            target, cr->rank_id, cr->comm, request);

  cr->n_elts[1] = 0;

  for (int i = 0; i < n_recv; i++) {

    unsigned char *r_ptr =   cr->buffer[0]
                           + _data_size(cr, cr->n_elts[0], cr->n_vals[0]);

    MPI_Irecv(r_ptr, recv_size[i*2 + 1], cr->mpi_type,
              target+i, target+i, cr->comm, request+1+i);

    cr->n_elts[0] += recv_size[i*2];
    cr->n_vals[0] += _comm_type_count_to_n_vals(cr,
                                                recv_size[i*2],
                                                recv_size[i*2 + 1]);

  }

  MPI_Waitall(n_recv+1, request, status);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers + 1, &t0, &t1);
}

/*----------------------------------------------------------------------------
 * Get destination id and associated data index for an indexed Crystal Router.
 *
 * If the destination id info is not provided in the Crystal Router,
 * it is assumed the given array is already valid, so only the index is
 * computed. If it is provided, the dest_id array is optional.
 *
 * parameters:
 *   <-- cr        pointer to associated Crystal Router
 *   <-> dest_id   pointer to destination id array, or NULL
 *   --> data_idx  pointer to data index
 *----------------------------------------------------------------------------*/

static void
_get_data_index_with_dest_id(cs_crystal_router_t  *cr,
                             cs_lnum_t             dest_id[],
                             cs_lnum_t             data_idx[])
{
  const size_t n_elts = cr->n_elts[0];

  const unsigned char *p_s = cr->buffer[0];

  if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {

    if (dest_id != NULL) {
      for (size_t i = 0; i < n_elts; i++) {
        const cs_lnum_t *cr_dest_id_p
          = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
        dest_id[i] = *cr_dest_id_p;
        const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
        data_idx[*cr_dest_id_p + 1] = *pn;
        p_s += cr->comp_size + cr->elt_size*(*pn);
      }
    }
    else {
      for (size_t i = 0; i < n_elts; i++) {
        const cs_lnum_t *cr_dest_id_p
          = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
        const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
        data_idx[*cr_dest_id_p + 1] = *pn;
        p_s += cr->comp_size + cr->elt_size*(*pn);
      }
    }

  }
  else {

    for (size_t i = 0; i < n_elts; i++) {
      const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
      data_idx[dest_id[i] + 1] = *pn;
      p_s += cr->comp_size + cr->elt_size*(*pn);
    }

  }

  /* Transform count to index */

  data_idx[0] = 0;
  for (size_t i = 0; i < n_elts; i++)
    data_idx[i+1] += data_idx[i];
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
 *   <-> dest_id    pointer to destination id array, or NULL
 *   --> src_rank   pointer to source rank array, or NULL
 *   --> src_id     pointer to source id array, or NULL
 *   --> data  pointer to (pointer to) destination data, or NULL
 *----------------------------------------------------------------------------*/

static void
_get_data_s_with_dest_id(cs_crystal_router_t   *cr,
                         cs_lnum_t             *dest_id,
                         int                   *src_rank,
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

  const unsigned char *cr_data_p = cr->buffer[0] + cr->elt_shift;
  unsigned char *data_p = (unsigned char *)data;

  /* If we do not have a dest_id buffer, use a local one */

  if (_dest_id == NULL) {
    BFT_MALLOC(_dest_id, n_elts, cs_lnum_t);
    cs_assert(cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID);
  }

  for (size_t c_id = 0; c_id < n_chunks; c_id++) {

    size_t i0 = c_id*chunk_size;
    size_t i1 = (c_id+1)*chunk_size;
    if (i1 > n_elts)
      i1 = n_elts;

    /* Extract _dest_id first, to use in other loops */

    if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
      for (size_t i = i0; i < i1; i++) {
        const cs_lnum_t *cr_dest_id_p
          = (const cs_lnum_t *)(  cr->buffer[0] + i*cr->comp_size
                                + cr->dest_id_shift);
        _dest_id[i] = *cr_dest_id_p;
      }
    }

    /* Optional source rank */

    if (src_rank != NULL) {
      for (size_t i = i0; i < i1; i++) {
        const int *cr_src_rank_p
          = (const cs_lnum_t *)(  cr->buffer[0] + i*cr->comp_size
                                + sizeof(int));
        src_rank[_dest_id[i]] = *cr_src_rank_p;
      }
    }

    /* Optional source id */

    if (src_id != NULL) {
      for (size_t i = i0; i < i1; i++) {
        const cs_lnum_t *cr_src_id_p
          = (const cs_lnum_t *)(  cr->buffer[0] + i*cr->comp_size
                                + cr->src_id_shift);
        src_id[_dest_id[i]] = *cr_src_id_p;
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

/*----------------------------------------------------------------------------
 * Get data elements associated with an indexed Crystal Router.
 *
 * This variant assumes no destination id is provided, so data is
 * copied in Crystal Router buffer order.
 *
 * parameters:
 *   <-- cr        pointer to associated Crystal Router
 *   --> src_rank  pointer to source rank array, or NULL
 *   --> src_id    pointer to source id array, or NULL
 *   --> data_idx  pointer to data index, or NULL
 *   --> data      pointer to destination data, or NULL
 *----------------------------------------------------------------------------*/

static void
_get_data_i(cs_crystal_router_t   *cr,
            int                   *src_rank,
            cs_lnum_t             *src_id,
            cs_lnum_t             *data_idx,
            void                  *data)
{
  const size_t n_elts = cr->n_elts[0];

  if (data_idx != NULL)
    data_idx[0] = 0;

  cs_lnum_t      s_idx = 0;
  unsigned char *data_p = (unsigned char *)data;
  unsigned const char *p_s = cr->buffer[0];

  /* Extract data */

  for (size_t i = 0; i < n_elts; i++) {

    /* Optional source rank */

    if (src_rank != NULL) {
      const int *cr_src_rank = (const int *)(p_s + sizeof(int));
      src_rank[i] = *cr_src_rank;
    }

    /* Optional source id */

    if (src_id != NULL) {
      const cs_lnum_t *cr_src_id = (const cs_lnum_t *)(p_s + cr->src_id_shift);
      src_id[i] = *cr_src_id;
    }

    /* Number of elements and optional data index */

    const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
    const cs_lnum_t n_sub = *pn;
    const size_t _sub_size = cr->elt_size*n_sub;

    if (data_idx != NULL)
      data_idx[i+1] = data_idx[i] + n_sub;

    /* Data */

    if (data_p != NULL) {
      unsigned const char *pe = p_s + cr->elt_shift;
      unsigned char *_data_p = data_p + s_idx*cr->elt_size;
      for (size_t j = 0; j < _sub_size; j++)
        _data_p[j] = pe[j];
    }

    s_idx += n_sub;
    p_s += cr->comp_size + _sub_size;

  }
}

/*----------------------------------------------------------------------------
 * Get data elements associated with an indexed Crystal Router.
 *
 * This variant assumes a destination id is provided, along with a matching
 * index. The destination id may either be provided or extracted from the
 * data if present; but the data index should always be provided (i.e.
 * precomputed), as it is required for indirect writes. If the destination
 * ids are both provided and present in the Crystal Router, they are assumed
 * to agree.
 *
 * parameters:
 *   <-- cr        pointer to associated Crystal Router
 *   <-- dest_id   pointer to destination id array, or NULL
 *   <-- data_idx  pointer to data index
 *   --> src_rank  pointer to source rank array, or NULL
 *   --> src_id    pointer to source id array, or NULL
 *   --> data      pointer to destination data, or NULL
 *----------------------------------------------------------------------------*/

static void
_get_data_i_with_dest_id(cs_crystal_router_t   *cr,
                         const cs_lnum_t        dest_id[],
                         const cs_lnum_t        data_idx[],
                         int                   *src_rank,
                         cs_lnum_t             *src_id,
                         void                  *data)
{
  const size_t n_elts = cr->n_elts[0];

  unsigned char *data_p = (unsigned char *)data;
  unsigned const char *p_s = cr->buffer[0];

  /* Extract data */

  for (size_t i = 0; i < n_elts; i++) {

    cs_lnum_t id;
    if (dest_id != NULL)
      id = dest_id[i];
    else {
      const cs_lnum_t *cr_dest_id
        = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
      id = *cr_dest_id;
    }

    /* Optional source rank */
    if (src_rank != NULL) {
      const int *cr_src_rank = (const int *)(p_s + sizeof(int));
      src_rank[id] = *cr_src_rank;
    }

    /* Optional source id */
    if (src_id != NULL) {
      const cs_lnum_t *cr_src_id = (const cs_lnum_t *)(p_s + cr->src_id_shift);
      src_id[id] = *cr_src_id;
    }

    /* Number of elements and optional data index */
    const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
    const cs_lnum_t n_sub = *pn;
    const size_t _sub_size = cr->elt_size*n_sub;
    assert(n_sub == data_idx[id+1] - data_idx[id]);

    /* Data */
    if (data_p != NULL) {
      unsigned const char *pe = p_s + cr->elt_shift;
      unsigned char *_data_p = data_p + data_idx[id]*cr->elt_size;
      for (size_t j = 0; j < _sub_size; j++)
        _data_p[j] = pe[j];
    }

    p_s += cr->comp_size + _sub_size;

  }
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

    unsigned char *p_s = cr->buffer[0] + i*cr->comp_size;
    unsigned char *pe = p_s + cr->elt_shift;
    unsigned const char *_psrc = _elt + i*cr->elt_size;

    int *pr = (int *)p_s;
    pr[0] = dest_rank[i];
    if (add_src_rank) {
      pr[1] = cr->rank_id;
    }

    if (add_dest_id) {
      cs_lnum_t *_cr_dest_id = (cs_lnum_t *)(p_s + cr->dest_id_shift);
      *_cr_dest_id = dest_id[i];
    }

    if (add_src_id) {
      cs_lnum_t *_cr_src_id = (cs_lnum_t *)(p_s + cr->src_id_shift);
      *_cr_src_id = (cs_lnum_t)i;
    }

    for (size_t j = 0; j < cr->elt_size; j++)
      pe[j] = _psrc[j];
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);

#if _CR_DEBUG_DUMP
  _dump(cr, "Crystal Router after creation.");
#endif

  return cr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a Crystal Router for indexed data.
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
 * \param[in]  elt_idx           element values start and past-the-last index
 * \param[in]  elt               element values
 * \param[in]  dest_id           element destination id, or NULL
 * \param[in]  dest_rank         destination rank for each element
 * \param[in]  comm              associated MPI communicator
 *
 * \return  pointer to new Crystal Router structure.
 */
/*----------------------------------------------------------------------------*/

cs_crystal_router_t *
cs_crystal_router_create_i(size_t            n_elts,
                           cs_datatype_t     datatype,
                           int               flags,
                           const cs_lnum_t  *elt_idx,
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

  cs_crystal_router_t *cr = _crystal_create_meta_i(n_elts,
                                                   datatype,
                                                   elt_idx,
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

    unsigned char *p_s =   cr->buffer[0] + i*cr->comp_size
                         + elt_idx[i]*cr->elt_size;
    int *pr = (int *)p_s;
    cs_lnum_t *pn = (cs_lnum_t *)(p_s + cr->n_vals_shift);
    unsigned char *pe = p_s + cr->elt_shift;
    unsigned const char *_psrc = _elt + elt_idx[i]*cr->elt_size;

    pr[0] = dest_rank[i];

    if (add_src_rank)
      pr[1] = cr->rank_id;

    if (add_dest_id) {
      unsigned char *pi = p_s + cr->dest_id_shift;
      unsigned const char *_p_dest_id = (const unsigned char *)(dest_id + i);
      for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
        pi[j] = _p_dest_id[j];
    }

    if (add_src_id) {
      const cs_lnum_t src_id = i;
      const unsigned char *_src_id = (const unsigned char *)(&src_id);
      unsigned char *pi = p_s + cr->src_id_shift;
      for (size_t j = 0; j < sizeof(cs_lnum_t); j++)
        pi[j] = _src_id[j];
    }

    cs_lnum_t n_sub = elt_idx[i+1] - elt_idx[i];
    pn[0] = n_sub;

    size_t sub_size = n_sub * cr->elt_size;
    for (size_t j = 0; j < sub_size; j++)
      pe[j] = _psrc[j];
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(_cr_timers, &t0, &t1);

#if _CR_DEBUG_DUMP
  _dump(cr, "Crystal Router after creation.");
#endif

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
      if (_cr->mpi_type != MPI_BYTE)
        MPI_Type_free(&(_cr->mpi_type));
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

    if (cr->n_vals_shift == 0)
      _crystal_partition_strided(cr, send_part, b_high);
    else
      _crystal_partition_indexed(cr, send_part, b_high);

#if _CR_DEBUG_DUMP
    {
      char comment[80];
      snprintf(comment, 79, "Crystal Router after partition, n_sub = %d",
               n_sub_ranks);
      _dump(cr, comment);
    }
#endif

    /* Send message to target process */

    _crystal_sendrecv(cr, target, n_recv);

    /* Ready for next exchange */

    if (cr->rank_id < b_high)
      n_sub_ranks = n_low;
    else {
      n_sub_ranks -= n_low;
      b_low = b_high;
    }

#if _CR_DEBUG_DUMP
    {
      char comment[80];
      if (n_sub_ranks > 1) {
        snprintf(comment, 79, "Crystal Router after sendrecv, n_sub = %d",
                 n_sub_ranks);
      }
      else {
        snprintf(comment, 79, "Crystal Router after exchange");
      }
      _dump(cr, comment);
    }
#endif

  }

  cr->n_elts[1] = 0;
  cr->buffer_size[1] = 0;
  BFT_FREE(cr->buffer[1]);

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

  if (cr != NULL) {

    if (cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {

      const size_t n_elts = cr->n_elts[0];
      cs_lnum_t dest_id_max = -1;

      if (cr->n_vals_shift == 0) {
        for (size_t i = 0; i < n_elts; i++) {
          unsigned const char *p_s = cr->buffer[0] + i*cr->comp_size;
          const cs_lnum_t *cr_dest_id
            = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
          if (cr_dest_id[0] > dest_id_max)
            dest_id_max = cr_dest_id[0];
        }
      }

      else {
        unsigned const char *p_s = cr->buffer[0];
        for (size_t i = 0; i < n_elts; i++) {
          const cs_lnum_t *cr_dest_id
            = (const cs_lnum_t *)(p_s + cr->dest_id_shift);
          if (cr_dest_id[0] > dest_id_max)
            dest_id_max = cr_dest_id[0];
          const cs_lnum_t *pn = (const cs_lnum_t *)(p_s + cr->n_vals_shift);
          p_s += cr->comp_size + cr->elt_size*pn[0];
        }
      }

      retval = dest_id_max + 1;

    }
    else
      retval = cr->n_elts[0];

  }

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
 * \param[in]       cr          pointer to associated Crystal Router
 * \param[out]      src_rank    pointer to (pointer to) source rank array,
 *                              or NULL
 * \param[in, out]  dest_id     pointer to (pointer to) destination id array,
 *                              or NULL
 * \param[out]      src_id      pointer to (pointer to) source id array, or NULL
 * \param[out]      data_index  pointer to (pointer to) destination index,
 *                              or NULL
 * \param[out]      data        pointer to (pointer to) destination data,
 *                              or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_crystal_router_get_data(cs_crystal_router_t   *cr,
                           int                  **src_rank,
                           cs_lnum_t            **dest_id,
                           cs_lnum_t            **src_id,
                           cs_lnum_t            **data_index,
                           void                 **data)
{
  cs_timer_t t0 = cs_timer_time();

  int *_src_rank = NULL;
  cs_lnum_t *_dest_id = NULL;
  cs_lnum_t *_src_id = NULL;
  cs_lnum_t *_data_index = NULL;
  unsigned char *_data = NULL;

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

  if (data_index != NULL && cr->n_vals_shift > 0) {
    _data_index = *data_index;
    if (_data_index == NULL) {
      BFT_MALLOC(_data_index, n_elts + 1, cs_lnum_t);
      *data_index = _data_index;
    }
  }

  if (data != NULL) {
    _data = *data;
    if (_data == NULL) {
      size_t data_size;
      if (cr->stride > 0)
        data_size = n_elts*cr->elt_size;
      else
        data_size = cr->n_vals[0]*cr->elt_size;
      BFT_MALLOC(_data, data_size, unsigned char);
      *data = _data;
    }
  }

  if (data != NULL && cr->stride > 0) {
    _data = *data;
    if (_data == NULL) {
      BFT_MALLOC(_data, n_elts*cr->elt_size, unsigned char);
      *data = _data;
    }
  }

  /* Now extract data */

  if (_dest_id != NULL || cr->flags & CS_CRYSTAL_ROUTER_USE_DEST_ID) {
    if (cr->n_vals_shift == 0)
      _get_data_s_with_dest_id(cr,
                               _dest_id,
                               _src_rank,
                               _src_id,
                               _data);
    else {
      if (cr->n_vals_shift > 0 && _data_index == NULL)
        BFT_MALLOC(_data_index, n_elts + 1, cs_lnum_t);

      _get_data_index_with_dest_id(cr, _dest_id, _data_index);

      _get_data_i_with_dest_id(cr,
                               _dest_id,
                               _data_index,
                               _src_rank,
                               _src_id,
                               _data);
    }
  }

  else {

    if (cr->n_vals_shift == 0)
      _get_data_s(cr, _src_rank, _src_id, _data);
    else
      _get_data_i(cr, _src_rank, _src_id, _data_index, _data);

  }

  if (dest_id == NULL && _dest_id != NULL)
    BFT_FREE(_dest_id);
  if (data_index == NULL && _data_index != NULL)
    BFT_FREE(_data_index);

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
