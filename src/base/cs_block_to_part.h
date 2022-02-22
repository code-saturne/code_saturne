#ifndef __CS_BLOCK_TO_PART_H__
#define __CS_BLOCK_TO_PART_H__

/*============================================================================
 * Convert between block distribution and general domain partition.
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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize block to partition distributor with block data using
 *        strided adjacency array.
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
 * The returned all to all distributor should be used in reverse
 * mode to copy data from the block to partition distribution.
 *
 * If the part_gnum array is requested (i.e. passed an non-NULL pointer),
 * the caller is responsible for freeing it.
 *
 * \param[in]   comm               communicator
 * \param[in]   block              block size and range info
 * \param[in]   adjacent_block     block info for adjacent entities
 * \param[in]   stride             stride of adjacency array (1 or 2)
 * \param[in]   adjacency          entity adjacency (1 to n numbering)
 * \param[in]   adjacent_ent_rank  destination rank for adjacent entities, or
 *                                 NULL if based on block size and range only.
 * \param[in]   default_rank       default rank in case there is no adjacency,
 *                                 or NULL
 * \param[out]  n_part_elts        number of elements in partition, or NULL
 * \param[out]  part_gnum          global elements numbers in partition,
 *                                 or NULL
 *
 * \return  initialized all to all block to partition distributor
 */
/*----------------------------------------------------------------------------*/

cs_all_to_all_t *
cs_block_to_part_create_by_adj_s(MPI_Comm               comm,
                                 cs_block_dist_info_t   block,
                                 cs_block_dist_info_t   adjacent_block,
                                 int                    stride,
                                 const cs_gnum_t        adjacency[],
                                 const int              adjacent_ent_rank[],
                                 const int              default_rank[],
                                 cs_lnum_t             *n_part_elts,
                                 cs_gnum_t            **part_gnum);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine local references from references to global numbers.
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
 * \param[in]   n_ents                 number of entities
 * \param[in]   base                   base numbering (typically 0 or 1)
 * \param[in]   global_list_size       size of global entity list
 * \param[in]   global_list_is_sorted  true if global entity list is
 *                                     guaranteed to be sorted
 * \param[in]   global_list            global entity list
 * \param[in]   global_number          entity global numbers (size: n_ents)
 * \param[out]  local_number           entity local numbers (size: n_ents)
 */
/*----------------------------------------------------------------------------*/

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
