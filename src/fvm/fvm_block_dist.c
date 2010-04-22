/*============================================================================
 * Utility functions for block distribution.
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
#include "fvm_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_block_dist.h"

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

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution.
 *
 * arguments:
 *   rank_id        <-- id of local rank (ignored in serial mode)
 *   n_ranks        <-- number of associated ranks
 *   min_rank_step  <-- minimum rank step between blocks
 *   min_block_size <-- minimum number of entities per block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

fvm_block_dist_info_t
fvm_block_dist_compute_sizes(int         rank_id,
                             int         n_ranks,
                             int         min_rank_step,
                             fvm_lnum_t  min_block_size,
                             fvm_gnum_t  n_g_ents)
{
  int _rank_id = rank_id;
  fvm_gnum_t _min_block_size = 1;
  fvm_gnum_t _block_size = 0;
  fvm_gnum_t _n_ranks = n_ranks;

  fvm_block_dist_info_t bi;

  /* Special case: only 1 rank */

  if (n_ranks == 1) {

    bi.gnum_range[0] = 1;
    bi.gnum_range[1] = n_g_ents + 1;
    bi.n_ranks = 1;
    bi.rank_step = 1;
    bi.block_size = n_g_ents;

    return bi;
  }

  /* Determine rank stepping if necessary */

  assert(rank_id > -1);

  bi.rank_step = 1;

  if (min_block_size > 1)
    _min_block_size = min_block_size;

  while (   n_g_ents/_n_ranks < _min_block_size
         && _n_ranks > 1
         && bi.rank_step < n_ranks) {
    bi.rank_step *= 2;
    _n_ranks = n_ranks / bi.rank_step;
  }
  if (bi.rank_step < min_rank_step)
    bi.rank_step = min_rank_step;
  if (bi.rank_step > n_ranks) {
    bi.rank_step = n_ranks;
    _n_ranks = 1;
  }

  if (rank_id % bi.rank_step == 0)
    _rank_id = rank_id/bi.rank_step;          /* non-empty block */
  else
    _rank_id = - (rank_id/bi.rank_step + 1);  /* empty block on this rank */

  /* Now determine block size and local range */

  _block_size = n_g_ents / _n_ranks;

  if (n_g_ents % _n_ranks)
    _block_size += 1;

  if (_rank_id > -1) {
    int i;
    for (i = 0; i < 2; i++) {
      fvm_gnum_t _g_rank = _rank_id + i;
      bi.gnum_range[i] = _g_rank*_block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }
  else {
    int i;
    fvm_gnum_t _g_rank = -_rank_id;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = _g_rank*_block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }

  bi.n_ranks = _n_ranks;
  bi.block_size = _block_size;

  return bi;
}

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

fvm_block_dist_info_t
fvm_block_dist_compute_sizes_nr(int         rank_id,
                                int         n_ranks,
                                int         n_block_ranks,
                                fvm_gnum_t  n_g_ents)
{
  int _rank_id = rank_id;

  fvm_gnum_t _block_size = 0;
  fvm_gnum_t _n_ranks = n_ranks;

  fvm_block_dist_info_t bi;

  /* Special case: only 1 rank */

  if (n_ranks == 1) {

    bi.gnum_range[0] = 1;
    bi.gnum_range[1] = n_g_ents + 1;
    bi.n_ranks = 1;
    bi.rank_step = 1;
    bi.block_size = n_g_ents;

    return bi;
  }

  /* Determine rank stepping if necessary */

  assert(rank_id > -1);

  _block_size = 0;
  _n_ranks = n_block_ranks;
  bi.rank_step = n_ranks / n_block_ranks;

  if (n_block_ranks < 1 || bi.rank_step > n_ranks) {
    bi.rank_step = n_ranks;
    _n_ranks = 1;
  }
  else if (bi.rank_step < 1) {
    bi.rank_step = 1;
    _n_ranks = n_ranks;
  }

  if (rank_id % bi.rank_step == 0)
    _rank_id = rank_id/bi.rank_step;          /* non-empty block */
  else
    _rank_id = - (rank_id/bi.rank_step + 1);  /* empty block on this rank */

  /* Now determine block size and local range */

  _block_size = n_g_ents / _n_ranks;

  if (n_g_ents % _n_ranks)
    _block_size += 1;

  if (_rank_id > -1) {
    int i;
    for (i = 0; i < 2; i++) {
      fvm_gnum_t _g_rank = _rank_id + i;
      bi.gnum_range[i] = _g_rank*_block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }
  else {
    int i;
    fvm_gnum_t _g_rank = -_rank_id;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = _g_rank*_block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }

  bi.n_ranks = _n_ranks;
  bi.block_size = _block_size;

  return bi;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

