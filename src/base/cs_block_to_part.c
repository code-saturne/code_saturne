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
#include "cs_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_block_to_part.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
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
                                 cs_gnum_t            **part_gnum)
{
  const cs_lnum_t n_ents = block.gnum_range[1] - block.gnum_range[0];

  int rank = -1;
  MPI_Comm_rank(comm, &rank);

  /* Determine an entity's adjacent entities, and use their
     global number to determine by which ranks they were read;
     using this knowledge, we can query these ranks to know to
     which ranks these adjacent entities were sent
     --------------------------------------------------------- */

  /* Count values to send and receive to each processor */

  cs_lnum_t _n_ents = n_ents*stride;

  int *query_rank;
  BFT_MALLOC(query_rank, _n_ents, int);

  for (cs_lnum_t j = 0; j < _n_ents; j++) {
    cs_gnum_t adj_g_num = adjacency[j];
    if (adj_g_num > 0) {
      int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                         * adjacent_block.rank_step;
      query_rank[j] = adj_ent_rank;
    }
    else
      query_rank[j] = rank; /* leave on current rank */
  }

  cs_all_to_all_t *qd = cs_all_to_all_create(_n_ents,
                                             0, /* flags */
                                             NULL,
                                             query_rank,
                                             comm);

  cs_all_to_all_transfer_dest_rank(qd, &query_rank);

  cs_gnum_t *adj_query
    = cs_all_to_all_copy_array(qd,
                               CS_GNUM_TYPE,
                               1,
                               false, /* reverse */
                               adjacency,
                               NULL);

  cs_lnum_t n_elts_query = cs_all_to_all_n_elts_dest(qd);

  int *sent_rank;
  BFT_MALLOC(sent_rank, n_elts_query, int);

  if (adjacent_ent_rank != NULL) {
    for (cs_lnum_t j = 0; j < n_elts_query; j++) {
      if (adj_query[j] > 0) {
        cs_lnum_t adj_l_id = (adj_query[j] - 1) % adjacent_block.block_size;
        sent_rank[j] = adjacent_ent_rank[adj_l_id];
      }
      else
        sent_rank[j] = -1;
    }
  }
  else {
    for (cs_lnum_t j = 0; j < n_elts_query; j++) {
      if (adj_query[j] > 0)
        sent_rank[j] = rank;
      else
        sent_rank[j] = -1;
    }
  }

  BFT_FREE(adj_query);

  /* Return communication */

  int *dest_rank
    = cs_all_to_all_copy_array(qd,
                               CS_INT_TYPE,
                               1,
                               true, /* reverse */
                               sent_rank,
                               NULL);

  BFT_FREE(sent_rank);

  cs_all_to_all_destroy(&qd);

  /* We now need to extract the ranks to which each entity will be sent,
     based on where its adjacent entities were sent (and thus where
     it will be needed).
     -------------------------------------------------------------------- */

  int  *send_rank;
  BFT_MALLOC(send_rank, _n_ents, int);
  cs_gnum_t *send_gnum;
  BFT_MALLOC(send_gnum, _n_ents, cs_gnum_t);

  cs_lnum_t n_send = 0;

  if (stride == 1) {
    for (cs_lnum_t j = 0; j < n_ents; j++) {
      int _send_rank = dest_rank[j];
      if (_send_rank != -1) {
        send_rank[n_send] = _send_rank;
        send_gnum[n_send] = block.gnum_range[0] + j;
        n_send++;
      }
      else if (default_rank != NULL) {
        send_rank[n_send] = default_rank[j];
        send_gnum[n_send] = block.gnum_range[0] + j;
        n_send++;
      }
    }
  }

  else if (stride == 2) {
    for (cs_lnum_t j = 0; j < n_ents; j++) {
      int _send_rank = -1;
      int prev_rank = -1;
      for (cs_lnum_t k = 0; k < 2; k++) {
        _send_rank = dest_rank[j*2 + k];
        if (_send_rank != -1 && _send_rank != prev_rank) {
          send_rank[n_send] = _send_rank;
          send_gnum[n_send] = block.gnum_range[0] + j;
          n_send++;
          prev_rank = _send_rank;
        }
      }
      if (prev_rank == -1 && default_rank != NULL) {
        send_rank[n_send] = default_rank[j];
        send_gnum[n_send] = block.gnum_range[0] + j;
        n_send++;
      }
    }
  }

  else
    bft_error(__FILE__, __LINE__, 0,
              "%s currently only allows stride 1 or 2, not %d.",
              __func__, (int)stride);

  BFT_FREE(dest_rank);

  /* Exchange send count and global number */

  /* Remark: we could save a global cs_all_to_all exchange by adding a
     source element list argument to the cs_all_to_all_create_*
     functions (such as the src_id in cs_crystal_router_create*),
     or adding a creation function with such a list;
     such a list is used internally for reverse cs_all_to_all
     exchanges, but allowing it would probably imply
     CS_ALL_TO_ALL_NO_REVERSE to avoid a too complex logic. */

  cs_all_to_all_t *d
    = cs_all_to_all_create(n_send,
                           CS_ALL_TO_ALL_ORDER_BY_SRC_RANK,
                           NULL,
                           send_rank,
                           comm);

  cs_gnum_t *recv_gnum
    = cs_all_to_all_copy_array(d,
                               CS_GNUM_TYPE,
                               1,
                               false, /* reverse */
                               send_gnum,
                               NULL);

  cs_lnum_t  _n_part_elts = cs_all_to_all_n_elts_dest(d);

  BFT_FREE(send_rank);
  BFT_FREE(send_gnum);

  cs_all_to_all_destroy(&d);

  /* Now build final distributor */

  d = cs_all_to_all_create_from_block(_n_part_elts,
                                      CS_ALL_TO_ALL_USE_DEST_ID,
                                      recv_gnum,
                                      block,
                                      comm);

  if (n_part_elts != NULL)
    *n_part_elts = _n_part_elts;

  if (part_gnum != NULL)
    *part_gnum = recv_gnum;
  else
    BFT_FREE(recv_gnum);

  /* Return initialized structure */

  return d;
}

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
                                 cs_lnum_t        local_number[])
{
  cs_lnum_t i;
  cs_lnum_t *order = NULL;
  cs_gnum_t *_g_list = NULL;
  const cs_gnum_t *g_list = global_list;

  if (n_ents == 0)
    return;

 #if defined(DEBUG) && !defined(NDEBUG)
  if (global_list_is_sorted) {
    for (i = 1; i < global_list_size; i++)
      assert(global_list[i] > global_list[i-1]);
  }
#endif

  if (global_list_is_sorted == false) {
    BFT_MALLOC(_g_list, global_list_size, cs_gnum_t);
    order = cs_order_gnum(NULL, global_list, global_list_size);
    for (i = 0; i < global_list_size; i++)
      _g_list[i] = global_list[order[i]];
    g_list = _g_list;
  }

  for (i = 0; i < n_ents; i++) {

    cs_lnum_t start_id = 0;
    cs_lnum_t end_id = global_list_size;

    const cs_gnum_t num_1 = global_number[i];

    /* Use binary search */

    while (start_id < end_id) {
      cs_lnum_t mid_id = start_id + ((end_id - start_id) / 2);
      if (g_list[mid_id] < num_1)
        start_id = mid_id + 1;
      else
        end_id = mid_id;  /* Can't be end_id = mid_id -1;
                             g_list[mid_id] >= num_1, so
                             end_id must not be < mid_id in case
                             g_list[mid_id] == num_1 */
    }

    /* start_id == end_id at this stage; */

    if (start_id < global_list_size && g_list[start_id] == num_1)
      local_number[i] = start_id + base;
    else
      local_number[i] = base - 1;

  }

  BFT_FREE(_g_list);

  if (order != NULL) {
    for (i = 0 ; i < n_ents ; i++)
      local_number[i] = order[local_number[i] - base] + base;
    BFT_FREE(order);
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
