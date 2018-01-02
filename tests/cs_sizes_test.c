/*============================================================================
 * Unit test for some FVM global numbering issues;
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_printf.h>

#include "cs_block_dist.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  int i, j;
  int rank_id[] = {0, 1024, 2048, 4095};
  int n_ranks = 4096;
  cs_gnum_t n_g_ents[] = {500, 1100000000, 2200000000,5400000000};

  cs_block_dist_info_t bi;

  for (i = 0; i < 4; i++) {

    int64_t g_ent_id = 0;
    int l_rank_id = 0;

    bft_printf("\ntest for n_g_ents = %llu\n",
               (unsigned long long)(n_g_ents[i]));
    for (j = 0; j < 4; j++) {
      bi = cs_block_dist_compute_sizes(rank_id[j],
                                       n_ranks,
                                       0,
                                       1024*1024,
                                       n_g_ents[i]);
      bft_printf("%d/%d [%llu %llu] block step %d, block size %d\n",
                 rank_id[j], n_ranks, (unsigned long long)(bi.gnum_range[0]),
                 (unsigned long long)(bi.gnum_range[1]),
                 (int)(bi.rank_step), (int)(bi.block_size));
    }

    g_ent_id = n_g_ents[i] - 1;
    l_rank_id = g_ent_id/bi.block_size * bi.rank_step;
    bft_printf("\nrank id for ent_id = %llu: %d\n",
               (unsigned long long)g_ent_id, l_rank_id);
  }

  exit(EXIT_SUCCESS);
}
