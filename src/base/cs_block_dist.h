#ifndef __CS_BLOCK_DIST_H__
#define __CS_BLOCK_DIST_H__

/*============================================================================
 * Definition of a block distribution.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Information structure for block size and entity range */

typedef struct {

  cs_gnum_t    gnum_range[2];  /* Start and past-the-end global numbers
                                  associated with local block */
  int          n_ranks;        /* Number of active ranks */
  int          rank_step;      /* Step between active block ranks
                                  (1 in basic case, > 1 if we seek to
                                  avoid too small buffers and agglomerate
                                  blocks on only a few ranks) */
  cs_lnum_t    block_size;     /* Basic block size */

} cs_block_dist_info_t;

/*=============================================================================
 * Public function prototypes
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

cs_block_dist_info_t
cs_block_dist_compute_sizes(int        rank_id,
                            int        n_ranks,
                            int        min_rank_step,
                            cs_lnum_t  min_block_size,
                            cs_gnum_t  n_g_ents);

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution
 * for a new global number of entities with a given number of active
 * ranks.
 *
 * arguments:
 *   rank_id        <-- id of local rank (ignored in serial mode)
 *   n_ranks        <-- number of associated ranks
 *   n_block_ranks  <-- number of ranks associated with a block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

cs_block_dist_info_t
cs_block_dist_compute_sizes_nr(int        rank_id,
                               int        n_ranks,
                               int        n_block_ranks,
                               cs_gnum_t  n_g_ents);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLOCK_DIST_H__ */
