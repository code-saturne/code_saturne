
/*============================================================================
 * Convert between general domain partition and block distribution.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2008-2009  EDF

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

/*----------------------------------------------------------------------------*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_config_defs.h"

#include "fvm_block_dist.h"
#include "fvm_order.h"
#include "fvm_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_part_to_block.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure used to redistribute data */

#if defined(HAVE_MPI)

struct _fvm_part_to_block_t {

  MPI_Comm     comm;            /* Associated MPI communicator */

  int          rank;            /* Local rank in communicator */
  int          n_ranks;         /* Number of ranks associated with
                                   communicator */

  fvm_part_to_block_info_t  bi; /* Associated block information */

  size_t       n_block_ents;    /* Number of entities to receive (this block) */
  size_t       n_part_ents;     /* Number of entities to send (partition) */
  size_t       recv_size;       /* Size of receive buffer for MPI_Alltoallv
                                   (send_size not necessary, as send_size
                                   should always be equal to n_part_ents,
                                   though elements may be assembled in a
                                   different order) */

  int         *send_count;      /* Send counts for MPI_Alltoall */
  int         *recv_count;      /* Receive counts for MPI_Alltoall */
  int         *send_displ;      /* Send displs for MPI_Alltoall */
  int         *recv_displ;      /* Receive displs for MPI_Alltoall */

  fvm_lnum_t  *recv_block_id;   /* Id in block of received entities */

  const fvm_gnum_t *global_ent_num; /* Shared global entity numbers */
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

static fvm_lnum_t
_compute_displ(int        n_ranks,
               const int  count[],
               int        displ[])
{
  int i;
  fvm_lnum_t total_count = 0;

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

static fvm_part_to_block_t *
_part_to_block_create(MPI_Comm comm)
{
  fvm_part_to_block_t *d;

  BFT_MALLOC(d, 1, fvm_part_to_block_t);

  d->comm = comm;

  MPI_Comm_rank(comm, &(d->rank));
  MPI_Comm_size(comm, &(d->n_ranks));

  memset(&(d->bi), 0, sizeof(fvm_part_to_block_info_t));

  d->n_block_ents = 0;
  d->n_part_ents = 0;
  d->recv_size = 0;

  d->send_count = NULL;
  d->recv_count = NULL;
  d->send_displ = NULL;
  d->recv_displ = NULL;

  d->recv_block_id = NULL;
  d->global_ent_num = NULL;

  return d;
}

/*----------------------------------------------------------------------------
 * Initialize partition to block distributor based on global element numbers,
 * using all to all communication.
 *
 * arguments:
 *   d    <-> partition to block distributor
 *   comm <-- communicator
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

static void
_init_alltoallv_by_gnum(fvm_part_to_block_t  *d,
                        MPI_Comm              comm)
{
  int i;
  size_t j;

  fvm_lnum_t *send_block_id = NULL;

  const int n_ranks = d->n_ranks;
  const int rank_step = d->bi.rank_step;
  const fvm_lnum_t block_size = d->bi.block_size;

  const fvm_gnum_t *global_ent_num = d->global_ent_num;

  /* Initialize send and receive counts */

  BFT_MALLOC(d->send_count, n_ranks, int);
  BFT_MALLOC(d->recv_count, n_ranks, int);
  BFT_MALLOC(d->send_displ, n_ranks, int);
  BFT_MALLOC(d->recv_displ, n_ranks, int);

  for (i = 0; i < d->n_ranks; i++)
    d->send_count[i] = 0;

  /* Count values to send and receive */

  for (j = 0; j < d->n_part_ents; j++) {
    int send_rank = ((global_ent_num[j] -1)/block_size)* rank_step;
    d->send_count[send_rank] += 1;
  }

  MPI_Alltoall(d->send_count, 1, MPI_INT, d->recv_count, 1, MPI_INT, comm);

  _compute_displ(n_ranks, d->send_count, d->send_displ);
  d->recv_size = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  /* Prepare list of local block ids of sent elements */

  BFT_MALLOC(d->recv_block_id, d->recv_size, fvm_lnum_t);

  BFT_MALLOC(send_block_id, d->n_part_ents, fvm_lnum_t);

  for (j = 0; j < d->n_part_ents; j++) {
    int64_t g_ent_id = global_ent_num[j] -1;
    int send_rank = (g_ent_id/block_size)*rank_step;
    send_block_id[d->send_displ[send_rank]] = g_ent_id % block_size;
    d->send_displ[send_rank] += 1;
  }

  /* Reset send_displ */

  for (i = 0; i < n_ranks; i++)
    d->send_displ[i] -= d->send_count[i];

  /* Exchange values */

  MPI_Alltoallv(send_block_id, d->send_count, d->send_displ, FVM_MPI_LNUM,
                d->recv_block_id, d->recv_count, d->recv_displ, FVM_MPI_LNUM,
                d->comm);

  BFT_FREE(send_block_id);
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
_init_gather_by_gnum(fvm_part_to_block_t  *d,
                     MPI_Comm              comm)
{
  size_t j;

  int send_count = d->n_part_ents;
  fvm_lnum_t *send_block_id = NULL;

  const int n_ranks = d->n_ranks;

  const fvm_gnum_t *global_ent_num = d->global_ent_num;

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
    BFT_MALLOC(d->recv_block_id, d->recv_size, fvm_lnum_t);

  BFT_MALLOC(send_block_id, d->n_part_ents, fvm_lnum_t);

  for (j = 0; j < d->n_part_ents; j++)
    send_block_id[j] = global_ent_num[j] -1;

  /* Exchange values */

  MPI_Gatherv(send_block_id, send_count, FVM_MPI_LNUM,
              d->recv_block_id, d->recv_count, d->recv_displ, FVM_MPI_LNUM,
              0, d->comm);

  BFT_FREE(send_block_id);
}

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition.
 *
 * arguments:
 *   d            <-- partition to block distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   part_values  <-- values in general domain partition
 *   block_values --> values in block distribution
 *----------------------------------------------------------------------------*/

static void
_copy_array_alltoallv(fvm_part_to_block_t   *d,
                      fvm_datatype_t         datatype,
                      int                    stride,
                      const void            *part_values,
                      void                  *block_values)
{
  int        i;
  size_t     j, k;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  size_t stride_size = fvm_datatype_size[datatype]*stride;
  MPI_Datatype mpi_type = fvm_datatype_to_mpi[datatype];

  unsigned char *_block_values = block_values;
  const unsigned char *_part_values = part_values;

  const int n_ranks = d->n_ranks;
  const int rank_step = d->bi.rank_step;
  const fvm_lnum_t base_block_size = d->bi.block_size;
  const size_t n_recv_ents = d->recv_size;

  const fvm_gnum_t *global_ent_num = d->global_ent_num;

  /* Adjust send and receive dimensions */

  if (stride > 1) {
    for (i = 0; i < n_ranks; i++) {
      d->send_count[i] *= stride;
      d->recv_count[i] *= stride;
      d->send_displ[i] *= stride;
      d->recv_displ[i] *= stride;
    }
  }

  /* Prepare MPI buffers */

  BFT_MALLOC(send_buf, d->n_part_ents*stride_size, unsigned char);

  /* Prepare list of element values to send */

  for (j = 0; j < d->n_part_ents; j++) {
    int64_t g_ent_id = global_ent_num[j] -1;
    int send_rank = (g_ent_id/base_block_size)*rank_step;
    size_t r_displ = j*stride_size;
    size_t w_displ = d->send_displ[send_rank]*stride_size;
    d->send_displ[send_rank] += 1;
    for (k = 0; k < stride_size; k++)
      send_buf[w_displ + k] = _part_values[r_displ + k];
  }

  /* Reset send_displ */

  for (i = 0; i < n_ranks; i++)
    d->send_displ[i] -= d->send_count[i];

  BFT_MALLOC(recv_buf, n_recv_ents*stride_size, unsigned char);

  /* Exchange values */

  MPI_Alltoallv(send_buf, d->send_count, d->send_displ, mpi_type,
                recv_buf, d->recv_count, d->recv_displ, mpi_type,
                d->comm);

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

  if (stride > 1) {
    for (i = 0; i < n_ranks; i++) {
      d->send_count[i] /= stride;
      d->recv_count[i] /= stride;
      d->send_displ[i] /= stride;
      d->recv_displ[i] /= stride;
    }
  }

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
_copy_array_gather(fvm_part_to_block_t   *d,
                   fvm_datatype_t         datatype,
                   int                    stride,
                   const void            *part_values,
                   void                  *block_values)
{
  int        i;
  size_t     j, k;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  int send_count = d->n_part_ents * stride;

  size_t stride_size = fvm_datatype_size[datatype]*stride;
  MPI_Datatype mpi_type = fvm_datatype_to_mpi[datatype];

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

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution.
 *
 * arguments:
 *   rank_id        <-- id of local rank
 *   n_ranks        <-- number of associated ranks
 *   min_rank_step  <-- minimum rank step between blocks
 *   min_block_size <-- minimum number of entities per block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

fvm_part_to_block_info_t
fvm_part_to_block_compute_sizes(int         rank_id,
                                int         n_ranks,
                                int         min_rank_step,
                                fvm_lnum_t  min_block_size,
                                fvm_gnum_t  n_g_ents)
{
  fvm_part_to_block_info_t bi;

  fvm_block_dist_info_t _bi= fvm_block_dist_compute_sizes(rank_id,
                                                          n_ranks,
                                                          min_rank_step,
                                                          min_block_size,
                                                          n_g_ents);

  bi.gnum_range[0] = _bi.gnum_range[0];
  bi.gnum_range[1] = _bi.gnum_range[1];
  bi.n_ranks = _bi.n_ranks;
  bi.rank_step = _bi.rank_step;
  bi.block_size = _bi.block_size;

  return bi;
}

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

fvm_part_to_block_t *
fvm_part_to_block_create_by_gnum(MPI_Comm                   comm,
                                 fvm_part_to_block_info_t   bi,
                                 fvm_lnum_t                 n_ents,
                                 const fvm_gnum_t           global_ent_num[])
{
  fvm_part_to_block_t *d = _part_to_block_create(comm);

  d->bi = bi;

  d->n_block_ents = bi.gnum_range[1] - bi.gnum_range[0];
  d->n_part_ents = n_ents;

  d->global_ent_num = global_ent_num;

  if (bi.n_ranks == 1)
    _init_gather_by_gnum(d, comm);
  else
    _init_alltoallv_by_gnum(d, comm);

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
fvm_part_to_block_destroy(fvm_part_to_block_t **d)
{
  fvm_part_to_block_t *_d = *d;

  BFT_FREE(_d->send_count);
  BFT_FREE(_d->recv_count);
  BFT_FREE(_d->send_displ);
  BFT_FREE(_d->recv_displ);

  BFT_FREE(_d->recv_block_id);

  BFT_FREE(*d);
}

/*----------------------------------------------------------------------------
 * Return number of entities associated with local partition
 *
 * arguments:
 *   d <-- distribtor helper
 *
 * returns:
 *   number of entities associated with distribution receive
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_part_to_block_get_n_part_ents(fvm_part_to_block_t *d)
{
  fvm_lnum_t retval = 0;

  if (d != NULL)
    retval = d->n_part_ents;

  return retval;
}

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition.
 *
 * arguments:
 *   d            <-- partition to block distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   part_values  <-- values in general domain partition
 *   block_values --> values in block distribution
 *----------------------------------------------------------------------------*/

void
fvm_part_to_block_copy_array(fvm_part_to_block_t   *d,
                             fvm_datatype_t         datatype,
                             int                    stride,
                             const void            *part_values,
                             void                  *block_values)
{
  if (d->bi.n_ranks == 1)
    _copy_array_gather(d,
                       datatype,
                       stride,
                       part_values,
                       block_values);
  else
    _copy_array_alltoallv(d,
                          datatype,
                          stride,
                          part_values,
                          block_values);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

