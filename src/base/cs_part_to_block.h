#ifndef __CS_PART_TO_BLOCK_H__
#define __CS_PART_TO_BLOCK_H__

/*============================================================================
 * Convert between general domain partition and block distribution.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_block_dist.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Opaque general domain partitioning to block distribution structure */

#if defined(HAVE_MPI)

typedef struct _cs_part_to_block_t  cs_part_to_block_t;

#endif

/*=============================================================================
 * Public function prototypes
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
cs_part_to_block_create_by_gnum(MPI_Comm              comm,
                                cs_block_dist_info_t  bi,
                                cs_lnum_t             n_ents,
                                const cs_gnum_t       global_ent_num[]);

/*----------------------------------------------------------------------------
 * Destroy a partition to block distributor structure.
 *
 * arguments:
 *   d <-> pointer to partition to block distributor structure pointer
 *----------------------------------------------------------------------------*/

void
cs_part_to_block_destroy(cs_part_to_block_t **d);

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
                               cs_gnum_t            global_ent_num[]);

/*----------------------------------------------------------------------------
 * Return number of entities associated with local partition
 *
 * arguments:
 *   d <-- distribtor helper
 *
 * returns:
 *   number of entities associated with distribution receive
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_part_to_block_get_n_part_ents(cs_part_to_block_t *d);

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
cs_part_to_block_copy_array(cs_part_to_block_t   *d,
                            cs_datatype_t         datatype,
                            int                   stride,
                            const void           *part_values,
                            void                 *block_values);

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
cs_part_to_block_copy_index(cs_part_to_block_t  *d,
                            const cs_lnum_t     *part_index,
                            cs_lnum_t           *block_index);

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
cs_part_to_block_copy_indexed(cs_part_to_block_t   *d,
                              cs_datatype_t         datatype,
                              const cs_lnum_t      *part_index,
                              const void           *part_val,
                              const cs_lnum_t      *block_index,
                              void                 *block_val);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PART_TO_BLOCK_H__ */
