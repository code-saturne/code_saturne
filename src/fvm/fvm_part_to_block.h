#ifndef __FVM_PART_TO_BLOCK_H__
#define __FVM_PART_TO_BLOCK_H__

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Information structure for block size and entity range */

typedef struct {

  fvm_gnum_t   gnum_range[2];  /* Start and past-the-end global numbers
                                  associated with local block */
  int          n_ranks;        /* Number of active ranks */
  int          rank_step;      /* Step between active block ranks
                                  (1 in basic case, > 1 if we seek to
                                  avoid too small buffers and agglomerate
                                  blocks on only a few ranks) */
  fvm_lnum_t   block_size;     /* Basic block size */

} fvm_part_to_block_info_t;

/* Opaque general domain partitioning to block distribution structure */

#if defined(HAVE_MPI)

typedef struct _fvm_part_to_block_t  fvm_part_to_block_t;

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a partitioned distribution.
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
                                fvm_gnum_t  n_g_ents);

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
                                 const fvm_gnum_t           global_ent_num[]);

/*----------------------------------------------------------------------------
 * Destroy a partition to block distributor structure.
 *
 * arguments:
 *   d <-> pointer to partition to block distributor structure pointer
 *----------------------------------------------------------------------------*/

void
fvm_part_to_block_destroy(fvm_part_to_block_t **d);

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
fvm_part_to_block_get_n_part_ents(fvm_part_to_block_t *d);

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
fvm_part_to_block_copy_array(fvm_part_to_block_t   *d,
                             fvm_datatype_t         datatype,
                             int                    stride,
                             const void            *part_values,
                             void                  *block_values);

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
 *                   (size: n_block_entities + 1)
 *----------------------------------------------------------------------------*/

void
fvm_part_to_block_copy_index(fvm_part_to_block_t  *d,
                             const fvm_lnum_t     *part_index,
                             fvm_lnum_t           *block_index);

/*----------------------------------------------------------------------------
 * Copy indexed data from general domain partition to block distribution.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- partition to block distributor
 *   datatype    <-- type of data considered
 *   part_index  <-- local index in general distribution
 *   part_val    <-- numbers in general  distribution
 *                   (size: part_index[n_part_ents])
 *   block_index --> local index in block distribution
 *   block_val   --> values in block distribution
 *                   (size: block_index[n_block_ents])
 *----------------------------------------------------------------------------*/

void
fvm_part_to_block_copy_indexed(fvm_part_to_block_t   *d,
                               fvm_datatype_t         datatype,
                               const fvm_lnum_t      *part_index,
                               const void            *part_val,
                               const fvm_lnum_t      *block_index,
                               void                  *block_val);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_PART_TO_BLOCK_H__ */
