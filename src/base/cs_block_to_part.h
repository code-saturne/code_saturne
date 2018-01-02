#ifndef __CS_BLOCK_TO_PART_H__
#define __CS_BLOCK_TO_PART_H__

/*============================================================================
 * Convert between block distribution and general domain partition.
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

/* Opaque block to general domain partitioning distribution structure */

#if defined(HAVE_MPI)

typedef struct _cs_block_to_part_t  cs_block_to_part_t;

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create block to partition distributor using entity destination rank array.
 *
 * arguments:
 *   comm     <-- communicator
 *   block    <-- block to partition range and size info
 *   ent_rank <-- destination rank for each entity
 *                (size: block.gnum_range[1] - block.gnum_range[0])
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

cs_block_to_part_t *
cs_block_to_part_create_by_rank(MPI_Comm              comm,
                                cs_block_dist_info_t  block,
                                int                   ent_rank[]);

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor with block data using
 * strided adjacency array.
 *
 * The adjacency array uses 1-n based global numbers. 0 values are
 * allowed and may be used to represent empty adjacencies.
 *
 * For example, in a face -> element adjacency relation, each face
 * is adjacent to 2 elements (thus a stride of 2), except for
 * boundary faces which are adjacent to only 1 element; in this case,
 * the adjacent element number for the exterior side of the face is 0.
 *
 * It is also possible to define a default destination rank,
 * so that elements with no adjacency are redistributed.
 * If the default rank for a given element is < 0, or no default
 * ranks are defined, elements with no adjacency are not distributed.
 *
 * arguments:
 *   comm              <-- communicator
 *   block             <-- block size and range info
 *   adjacent_block    <-- block info for adjacent entities
 *   stride            <-- stride of adjacency array
 *   adjacency         <-- entity adjacency (1 to n numbering)
 *   adjacent_ent_rank <-- destination rank for adjacent entities, or
 *                         NULL if based on block size and range only.
 *   default_rank      <-- default rank in case there is no adjacency, or NULL
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

cs_block_to_part_t *
cs_block_to_part_create_by_adj_s(MPI_Comm              comm,
                                 cs_block_dist_info_t  block,
                                 cs_block_dist_info_t  adjacent_block,
                                 int                   stride,
                                 cs_gnum_t             adjacency[],
                                 int                   adjacent_ent_rank[],
                                 int                   default_rank[]);

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor for entities adjacent to
 * already distributed entities.
 *
 * arguments:
 *   comm      <-- communicator
 *   bi        <-- block size and range info
 *   adj_bi    <-- block info for adjacent entities
 *   adjacency <-- entity adjacency (1 to n numbering)
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

cs_block_to_part_t *
cs_block_to_part_create_adj(MPI_Comm              comm,
                            cs_block_dist_info_t  adj_bi,
                            size_t                adjacency_size,
                            const cs_gnum_t       adjacency[]);

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor based global element numbers
 * for partitioned data.
 *
 * arguments:
 *   comm           <-- communicator
 *   bi             <-- block size and range info
 *   n_ents         <-- number of elements in partition
 *   global_ent_num <-- global entity numbers (in partition)
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

cs_block_to_part_t *
cs_block_to_part_create_by_gnum(MPI_Comm              comm,
                                cs_block_dist_info_t  bi,
                                cs_lnum_t             n_ents,
                                const cs_gnum_t       global_ent_num[]);

/*----------------------------------------------------------------------------
 * Destroy a block to partition distributor structure.
 *
 * arguments:
 *   d <-> pointer to block to partition distributor structure pointer
 *----------------------------------------------------------------------------*/

void
cs_block_to_part_destroy(cs_block_to_part_t  **d);

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
cs_block_to_part_get_n_part_ents(cs_block_to_part_t  *d);

/*----------------------------------------------------------------------------
 * Transfer a block to partition distributor's associated global numbering.
 *
 * The pointer to the global number array is returned, and ownership
 * of this array is given to the caller.
 *
 * arguments:
 *   d <-> pointer to block to partition distributor structure pointer
 *
 * returns:
 *   pointer to receiver global numbering, or NULL if the block to
 *   domain partition distributor was not the owner of this array.
 *----------------------------------------------------------------------------*/

cs_gnum_t *
cs_block_to_part_transfer_gnum(cs_block_to_part_t  *d);

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition.
 *
 * arguments:
 *   d            <-- block to partition distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   block_values --> values in block distribution
 *   part_values  --> values in general domain partition
 *----------------------------------------------------------------------------*/

void
cs_block_to_part_copy_array(cs_block_to_part_t  *d,
                            cs_datatype_t        datatype,
                            int                  stride,
                            const void          *block_values,
                            void                *part_values);

/*----------------------------------------------------------------------------
 * Copy a local index from block distribution to general domain partition.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- block to partition distributor
 *   block_index <-- local index in block distribution
 *   part_index  --> local index in general partition distribution
 *                   (size: n_part_entities + 1)
 *----------------------------------------------------------------------------*/

void
cs_block_to_part_copy_index(cs_block_to_part_t  *d,
                            const cs_lnum_t     *block_index,
                            cs_lnum_t           *part_index);

/*----------------------------------------------------------------------------
 * Copy indexed data from block distribution to general domain partition.
 *
 * arguments:
 *   d           <-- block to partition distributor
 *   datatype    <-- type of data considered
 *   block_index <-- local index in block distribution
 *   block_val   <-- values in block distribution
 *                   (size: block_index[n_block_ents])
 *   part_index  --> local index in general distribution
 *   part_val    --> numbers in general  distribution
 *                   (size: part_index[n_part_ents])
 *----------------------------------------------------------------------------*/

void
cs_block_to_part_copy_indexed(cs_block_to_part_t  *d,
                              cs_datatype_t        datatype,
                              const cs_lnum_t     *block_index,
                              const void          *block_val,
                              const cs_lnum_t     *part_index,
                              void                *part_val);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Determine local references from references to global numbers.
 *
 * This is based on finding the local id of a given global number
 * using a binary search.
 *
 * Global numbers use a 1 to n numbering, while local numbers use a
 * 0+base to n-1+base numbering. If an entity's global number does not
 * appear in the global list, base-1 is assigned for that entity's
 * local list.
 *
 * If list contains duplicate values, any local id having a multiple
 * global number (i.e not necessarily the smallest one) may be
 * assigned to the corresponding local_number[] entry.
 *
 * arguments:
 *   n_ents                <-- number of entities
 *   base                  <-- base numbering (typically 0 or 1)
 *   global_list_size      <-- size of global entity list
 *   global_list_is_sorted <-- true if global entity list is guaranteed
 *                             to be sorted
 *   global_list           <-- global entity list
 *   global_number         <-- entity global numbers
 *                             (size: n_ents)
 *   local_number          --> entity local numbers
 *                             (size: n_ents)
 *----------------------------------------------------------------------------*/

void
cs_block_to_part_global_to_local(cs_lnum_t        n_ents,
                                 cs_lnum_t        base,
                                 cs_lnum_t        global_list_size,
                                 bool             global_list_is_sorted,
                                 const cs_gnum_t  global_list[],
                                 const cs_gnum_t  global_number[],
                                 cs_lnum_t        local_number[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLOCK_TO_PART_H__ */
