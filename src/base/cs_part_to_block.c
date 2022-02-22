/*============================================================================
 * Convert between general domain partition and block distribution.
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
#include "bft_error.h"

#include "cs_all_to_all.h"
#include "cs_block_dist.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure used to redistribute data */

#if defined(HAVE_MPI)

struct _cs_part_to_block_t {

  MPI_Comm     comm;            /* Associated MPI communicator */

  int          rank;            /* Local rank in communicator */
  int          n_ranks;         /* Number of ranks associated with
                                   communicator */

  cs_block_dist_info_t  bi;     /* Associated block information */

  cs_all_to_all_t      *d;      /* Associated all-to-all distributor */

  size_t       n_block_ents;    /* Number of entities to receive (this block) */
  size_t       n_part_ents;     /* Number of entities to send (partition) */
  size_t       recv_size;       /* Size of receive buffer for MPI_Gatherv
                                   (send_size not necessary, as send_size
                                   should always be equal to n_part_ents,
                                   though elements may be assembled in a
                                   different order) */

  int         *recv_count;      /* Receive counts for MPI_Gatherv */
  int         *recv_displ;      /* Receive displs for MPI_Gatherv */

  int         *block_rank_id;   /* Block id for each part entity
                                   (NULL if based on global_ent_num) */
  cs_lnum_t   *send_block_id;   /* Id in block of sent entities */
  cs_lnum_t   *recv_block_id;   /* Id in block of received entities */

  const cs_gnum_t   *global_ent_num;  /* Shared global entity numbers */
  cs_gnum_t         *_global_ent_num; /* Private global entity numbers */
};

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
 *   count   <-- number of entities per rank (size: n_ranks)
 *   displ   --> entity displacement in cumulative array (size: n_ranks)
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
 * Create distribution helper structure.
 *
 * Send and receive counts and displacements are allocated, but not
 * fully initialized at this point: only the send count is set to zero.
 *
 * arguments:
 *   comm <-- communicator
 *
 * returns:
 *   empty communicator structure
 *----------------------------------------------------------------------------*/

static cs_part_to_block_t *
_part_to_block_create(MPI_Comm comm)
{
  cs_part_to_block_t *d;

  BFT_MALLOC(d, 1, cs_part_to_block_t);

  d->comm = comm;

  MPI_Comm_rank(comm, &(d->rank));
  MPI_Comm_size(comm, &(d->n_ranks));

  memset(&(d->bi), 0, sizeof(cs_block_dist_info_t));

  d->d = NULL;

  d->n_block_ents = 0;
  d->n_part_ents = 0;
  d->recv_size = 0;

  d->recv_count = NULL;
  d->recv_displ = NULL;

  d->block_rank_id = NULL;
  d->send_block_id = NULL;
  d->recv_block_id = NULL;
  d->global_ent_num = NULL;
  d->_global_ent_num = NULL;

  return d;
}

/*----------------------------------------------------------------------------
 * Initialize partition to block distributor based on global element numbers,
 * using gather to rank 0 when only one block is active.
 *
 * arguments:
 *   d    <-> partition to block distributor
 *   comm <-- communicator
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

static void
_init_gather_by_gnum(cs_part_to_block_t  *d,
                     MPI_Comm             comm)
{
  size_t j;

  int send_count = d->n_part_ents;
  cs_lnum_t *send_block_id = NULL;

  const int n_ranks = d->n_ranks;

  const cs_gnum_t *global_ent_num = d->global_ent_num;

  /* Initialize send and receive counts */

  if (d->rank == 0) {
    BFT_MALLOC(d->recv_count, n_ranks, int);
    BFT_MALLOC(d->recv_displ, n_ranks, int);
  }

  /* Count values to send and receive */

  MPI_Gather(&send_count, 1, MPI_INT, d->recv_count, 1, MPI_INT, 0, comm);

  if (d->rank == 0)
    d->recv_size = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  /* Prepare list of local block ids of sent elements */

  if (d->rank == 0)
    BFT_MALLOC(d->recv_block_id, d->recv_size, cs_lnum_t);

  BFT_MALLOC(send_block_id, d->n_part_ents, cs_lnum_t);

  for (j = 0; j < d->n_part_ents; j++)
    send_block_id[j] = global_ent_num[j] -1;

  /* Exchange values */

  MPI_Gatherv(send_block_id, send_count, CS_MPI_LNUM,
              d->recv_block_id, d->recv_count, d->recv_displ, CS_MPI_LNUM,
              0, d->comm);

  BFT_FREE(send_block_id);
}

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition,
 * using gather to rank 0 when only one block is active.
 *
 * arguments:
 *   d            <-- partition to block distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   part_values  <-- values in general domain partition
 *   block_values --> values in block distribution
 *----------------------------------------------------------------------------*/

static void
_copy_array_gatherv(cs_part_to_block_t  *d,
                    cs_datatype_t        datatype,
                    int                  stride,
                    const void          *part_values,
                    void                *block_values)
{
  int        i;
  size_t     j, k;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  int send_count = d->n_part_ents * stride;

  size_t stride_size = cs_datatype_size[datatype]*stride;
  MPI_Datatype mpi_type = cs_datatype_to_mpi[datatype];

  unsigned char *_block_values = block_values;

  const int n_ranks = d->n_ranks;
  const size_t n_recv_ents = d->recv_size;

  /* Adjust send and receive dimensions */

  if (stride > 1 && d->rank == 0) {
    for (i = 0; i < n_ranks; i++) {
      d->recv_count[i] *= stride;
      d->recv_displ[i] *= stride;
    }
  }

  BFT_MALLOC(recv_buf, n_recv_ents*stride_size, unsigned char);

  BFT_MALLOC(send_buf, d->n_part_ents*stride_size, unsigned char);
  if (d->n_part_ents > 0)
    memcpy(send_buf, part_values, d->n_part_ents*stride_size);

  /* Exchange values */

  MPI_Gatherv(send_buf, send_count, mpi_type,
              recv_buf, d->recv_count, d->recv_displ, mpi_type,
              0, d->comm);

  /* Distribute received values */

  for (j = 0; j < n_recv_ents; j++) {

    size_t r_displ = j*stride_size;
    size_t w_displ = d->recv_block_id[j]*stride_size;

    for (k = 0; k < stride_size; k++)
      _block_values[w_displ + k] = recv_buf[r_displ + k];
  }

  /* Cleanup */

  BFT_FREE(recv_buf);
  BFT_FREE(send_buf);

  /* Reset send and receive dimensions */

  if (stride > 1 && d->rank == 0) {
    for (i = 0; i < n_ranks; i++) {
      d->recv_count[i] /= stride;
      d->recv_displ[i] /= stride;
    }
  }
}

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition,
 * using gather to rank 0 when only one block is active.
 *
 * arguments:
 *   d           <-- partition to block distributor
 *   part_index  <-- local index in general partition distribution
 *                   (size: n_part_entities + 1)
 *   block_index --> local index in block distribution
 *                   (size: n_block_entities + 1)
 *----------------------------------------------------------------------------*/

static void
_copy_index_gatherv(cs_part_to_block_t  *d,
                    const cs_lnum_t     *part_index,
                    cs_lnum_t           *block_index)
{
  size_t     j;

  cs_lnum_t *send_buf = NULL;
  cs_lnum_t *recv_buf = NULL;

  int send_count = d->n_part_ents;

  const size_t n_recv_ents = d->recv_size;

  /* Prepare MPI buffers */

  BFT_MALLOC(send_buf, d->n_part_ents, cs_lnum_t);

  /* Prepare list of element values to send */

  for (j = 0; j < d->n_part_ents; j++)
    send_buf[j] = part_index[j+1] - part_index[j];

  BFT_MALLOC(recv_buf, n_recv_ents, cs_lnum_t);

  /* Exchange values */

  MPI_Gatherv(send_buf, send_count, CS_MPI_LNUM,
              recv_buf, d->recv_count, d->recv_displ, CS_MPI_LNUM,
              0, d->comm);

  /* Distribute received values */

  if (block_index != NULL) {

    for (j = 0; j < d->n_block_ents+1; j++)
      block_index[j] = 0;

    for (j = 0; j < n_recv_ents; j++) {
      assert(   block_index[d->recv_block_id[j]+1] == 0
             || block_index[d->recv_block_id[j]+1] == recv_buf[j]);
      block_index[d->recv_block_id[j]+1] = recv_buf[j];
    }

    /* Transform count to index */

    for (j = 0; j < d->n_block_ents; j++)
      block_index[j+1] += block_index[j];
  }

  /* Cleanup */

  BFT_FREE(recv_buf);
  BFT_FREE(send_buf);
}

/*----------------------------------------------------------------------------
 * Copy indexed data from general domain partition to block distribution,
 * using gather to rank 0 when only one block is active.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- partition to block distributor
 *   datatype    <-- type of data considered
 *   part_index  <-- local index in general distribution
 *   part_val    <-- numbers in general  distribution
 *                   (size: send_index[n_part_ents])
 *   block_index --> local index in block distribution
 *   block_val   --> values in block distribution
 *                   (size: recv_index[n_block_ents])
 *----------------------------------------------------------------------------*/

static void
_copy_indexed_gatherv(cs_part_to_block_t  *d,
                      cs_datatype_t        datatype,
                      const cs_lnum_t     *part_index,
                      const void          *part_val,
                      const cs_lnum_t     *block_index,
                      void                *block_val)
{
  int    i, l;
  size_t j, k;

  size_t w_displ = 0, r_displ = 0;
  size_t recv_size = 0;
  int  send_count = 0;

  int  *recv_count = NULL;
  int  *recv_displ = NULL;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  const unsigned char *_part_val = part_val;
  unsigned char *_block_val = block_val;

  size_t type_size = cs_datatype_size[datatype];
  MPI_Datatype mpi_type = cs_datatype_to_mpi[datatype];

  const int n_ranks = d->n_ranks;
  const size_t n_recv_ents = d->recv_size;

  /* Build send and receive counts */
  /*-------------------------------*/

  if (d->rank == 0) {
    BFT_MALLOC(recv_count, n_ranks, int);
    BFT_MALLOC(recv_displ, n_ranks, int);
    for (i = 0; i < n_ranks; i++)
      recv_count[i] = 0;
  }

  /* Prepare count of element values to send */

  for (j = 0; j < d->n_part_ents; j++)
    send_count += part_index[j+1] - part_index[j];

  /* Prepare count of element values to receive */

  if (d->rank == 0) {
    k = 0;
    for (i = 0; i < n_ranks; i++) {
      for (l = 0; l < d->recv_count[i]; l++) {
        w_displ = d->recv_block_id[k++];
        recv_count[i] += block_index[w_displ + 1] - block_index[w_displ];
      }
    }
    recv_size = _compute_displ(n_ranks, recv_count, recv_displ);
  }

  /* Build send and receive buffers */
  /*--------------------------------*/

  if (d->rank == 0)
    BFT_MALLOC(recv_buf, recv_size*type_size, unsigned char);

  BFT_MALLOC(send_buf, send_count * type_size, unsigned char);

  w_displ = 0;
  for (j = 0; j < d->n_part_ents; j++) {
    size_t ent_size = (part_index[j+1] - part_index[j]);
    r_displ = part_index[j]*type_size;
    for (k = 0; k < ent_size*type_size; k++)
      send_buf[w_displ + k] = _part_val[r_displ + k];
    w_displ += ent_size*type_size;
  }

  assert(w_displ == send_count*type_size);

  /* Exchange values */

  MPI_Gatherv(send_buf, send_count, mpi_type,
              recv_buf, recv_count, recv_displ, mpi_type,
              0, d->comm);

  BFT_FREE(send_buf);

  /* Distribute received values */

  if (block_index != NULL) {

    r_displ = 0;

    for (j = 0; j < n_recv_ents; j++) {

      size_t block_id = d->recv_block_id[j];
      size_t ent_size =   (block_index[block_id+1] - block_index[block_id])
                        * type_size;
      w_displ = block_index[block_id] * type_size;

      for (k = 0; k < ent_size; k++)
        _block_val[w_displ++] = recv_buf[r_displ++];
    }

    assert(r_displ == recv_size*type_size);
  }

  /* Cleanup */

  if (d->rank == 0) {
    BFT_FREE(recv_buf);
    BFT_FREE(recv_count);
    BFT_FREE(recv_displ);
  }
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize partition to block distributor based on global entity numbers.
 *
 * arguments:
 *   comm           <-- communicator
 *   bi             <-- block size and range info
 *   n_ents         <-- number of elements in partition
 *   global_ent_num <-- global entity numbers
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

cs_part_to_block_t *
cs_part_to_block_create_by_gnum(MPI_Comm               comm,
                                cs_block_dist_info_t   bi,
                                cs_lnum_t              n_ents,
                                const cs_gnum_t        global_ent_num[])
{
  cs_part_to_block_t *d = _part_to_block_create(comm);

  d->bi = bi;

  d->n_block_ents = bi.gnum_range[1] - bi.gnum_range[0];
  d->n_part_ents = n_ents;

  d->global_ent_num = global_ent_num;

  if (bi.n_ranks == 1)
    _init_gather_by_gnum(d, comm);
  else {
    int flags = CS_ALL_TO_ALL_USE_DEST_ID | CS_ALL_TO_ALL_NO_REVERSE;
    d->d = cs_all_to_all_create_from_block(n_ents,
                                           flags,
                                           global_ent_num,
                                           bi,
                                           comm);
  }

  /* Return initialized structure */

  return d;
}

/*----------------------------------------------------------------------------
 * Destroy a partition to block distributor structure.
 *
 * arguments:
 *   d <-> pointer to partition to block distributor structure pointer
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_destroy(cs_part_to_block_t  **d)
{
  cs_part_to_block_t *_d = *d;

  if (_d->d != NULL)
    cs_all_to_all_destroy(&(_d->d));

  BFT_FREE(_d->recv_count);
  BFT_FREE(_d->recv_displ);

  BFT_FREE(_d->block_rank_id);
  BFT_FREE(_d->send_block_id);
  BFT_FREE(_d->recv_block_id);

  if (_d->_global_ent_num != NULL)
    BFT_FREE(_d->_global_ent_num);

  BFT_FREE(*d);
}

/*----------------------------------------------------------------------------
 * Transfer ownership of global entity numbers to a block distributor.
 *
 * The global_ent_num[] array should be the same as the one used
 * for the creation of the block distributor.
 *
 * arguments:
 *   d              <-- distributor helper
 *   global_ent_num <-> global entity numbers
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_transfer_gnum(cs_part_to_block_t  *d,
                               cs_gnum_t            global_ent_num[])
{
  assert(d->global_ent_num == global_ent_num);

  d->_global_ent_num = global_ent_num;
}

/*----------------------------------------------------------------------------
 * Return number of entities associated with local partition
 *
 * arguments:
 *   d <-- distributor helper
 *
 * returns:
 *   number of entities associated with distribution receive
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_part_to_block_get_n_part_ents(cs_part_to_block_t  *d)
{
  cs_lnum_t retval = 0;

  if (d != NULL)
    retval = d->n_part_ents;

  return retval;
}

/*----------------------------------------------------------------------------
 * Copy array data from general domain partition to block distribution.
 *
 * arguments:
 *   d            <-- partition to block distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   part_values  <-- values in general domain partition
 *   block_values --> values in block distribution
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_copy_array(cs_part_to_block_t  *d,
                            cs_datatype_t        datatype,
                            int                  stride,
                            const void          *part_values,
                            void                *block_values)
{
  if (d->bi.n_ranks == 1)
    _copy_array_gatherv(d,
                        datatype,
                        stride,
                        part_values,
                        block_values);
  else
    cs_all_to_all_copy_array(d->d,
                             datatype,
                             stride,
                             false, /* reverse */
                             part_values,
                             block_values);
}

/*----------------------------------------------------------------------------
 * Copy local index from general domain partition to block distribution.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- partition to block distributor
 *   part_index  <-- local index in general partition distribution
 *                   (size: n_part_entities + 1)
 *   block_index --> local index in block distribution
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_copy_index(cs_part_to_block_t  *d,
                            const cs_lnum_t     *part_index,
                            cs_lnum_t           *block_index)
{
  if (d->bi.n_ranks == 1)
    _copy_index_gatherv(d, part_index, block_index);
  else
    cs_all_to_all_copy_index(d->d,
                             false, /* reverse */
                             part_index,
                             block_index);
}

/*----------------------------------------------------------------------------
 * Copy indexed data from general domain partition to block distribution,
 * using gather to rank 0 when only one block is active.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- partition to block distributor
 *   datatype    <-- type of data considered
 *   part_index  <-- local index in general distribution
 *   part_val    <-- numbers in general  distribution
 *                   (size: send_index[n_part_ents])
 *   block_index --> local index in block distribution
 *   block_val   --> values in block distribution
 *                   (size: recv_index[n_block_ents])
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_copy_indexed(cs_part_to_block_t  *d,
                              cs_datatype_t        datatype,
                              const cs_lnum_t     *part_index,
                              const void          *part_val,
                              const cs_lnum_t     *block_index,
                              void                *block_val)
{
  if (d->bi.n_ranks == 1)
    _copy_indexed_gatherv(d,
                          datatype,
                          part_index,
                          part_val,
                          block_index,
                          block_val);
  else
    cs_all_to_all_copy_indexed(d->d,
                               datatype,
                               false, /* reverse */
                               part_index,
                               part_val,
                               block_index,
                               block_val);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
