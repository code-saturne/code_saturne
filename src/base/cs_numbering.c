/*============================================================================
 * Numbering information for vectorization or multithreading
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_defs.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_numbering.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names for numbering types */

const char  *cs_numbering_type_name[] = {N_("default"),
                                         N_("vectorization"),
                                         N_("threads")};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a default numbering information structure.
 *
 * parameters:
 *   n_elts  <-- number of associated elements
 *
 * returns:
 *   pointer to created cs_numbering_t structure
 *---------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_default(cs_lnum_t  n_elts)
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_DEFAULT;

  numbering->vector_size = 1;

  numbering->n_threads = 1;
  numbering->n_groups = 1;

  numbering->n_no_adj_halo_groups = 0;
  numbering->n_no_adj_halo_elts = 0;

  BFT_MALLOC(numbering->group_index, 2, cs_lnum_t);
  numbering->group_index[0] = 0;
  numbering->group_index[1] = n_elts;

  return numbering;
}

/*----------------------------------------------------------------------------
 * Create a numbering information structure in case of vectorization.
 *
 * parameters:
 *   n_elts      <-- number of associated elements
 *   vector_size <-- vector size used for this vectorization
 *
 * returns:
 *   pointer to created cs_numbering_t structure
 *---------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_vectorized(cs_lnum_t  n_elts,
                               int        vector_size)
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_VECTORIZE;

  numbering->vector_size = vector_size;

  numbering->n_threads = 1;
  numbering->n_groups = 1;

  BFT_MALLOC(numbering->group_index, 2, cs_lnum_t);
  numbering->group_index[0] = 0;
  numbering->group_index[1] = n_elts;

  return numbering;
}

/*----------------------------------------------------------------------------
 * Create a numbering information structure in case of threading.
 *
 * parameters:
 *   n_threads   <-- number of threads
 *   n_groups    <-- number of groups
 *   group_index <-- group_index[thread_id*group_id*2 + group_id*2] and
 *                   group_index[thread_id*group_id*2 + group_id*2 +1] define
 *                   the start and end ids for entities in a given group and
 *                   thread; (size: n_groups *2 * n_threads)
 *
 * returns:
 *   pointer to created cs_numbering_t structure
 *---------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_threaded(int        n_threads,
                             int        n_groups,
                             cs_lnum_t  group_index[])
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_THREADS;

  numbering->vector_size = 1;

  numbering->n_threads = n_threads;
  numbering->n_groups = n_groups;

  BFT_MALLOC(numbering->group_index, n_threads*2*n_groups, cs_lnum_t);

  memcpy(numbering->group_index,
         group_index,
         (n_threads*2*n_groups) * sizeof(cs_lnum_t));

  return numbering;
}

/*----------------------------------------------------------------------------
 * Destroy a numbering information structure.
 *
 * parameters:
 *   numbering <-> pointer to cs_numbering_t structure pointer (or NULL)
 *---------------------------------------------------------------------------*/

void
cs_numbering_destroy(cs_numbering_t  **numbering)
{
  if (*numbering != NULL) {

    cs_numbering_t  *_n = *numbering;

    BFT_FREE(_n->group_index);

    BFT_FREE(*numbering);
  }
}

/*----------------------------------------------------------------------------
 * Dump a cs_numbering_t structure.
 *
 * parameters:
 *   numbering <-- pointer to cs_numbering_t structure (or NULL)
 *---------------------------------------------------------------------------*/

void
cs_numbering_dump(const cs_numbering_t  *numbering)
{
  int  i, j;

  if (numbering == NULL) {
    bft_printf("\n  Numbering: nil (default)\n");
    return;
  }

  bft_printf("\n  Numbering:           %p\n"
             "  type:                  %s\n"
             "  vector_size:           %d\n"
             "  n_threads:             %d\n"
             "  n_groups:              %d\n"
             "  n_no_adj_halo_groups:  %d\n"
             "  n_no_adj_halo_elts:    %ld\n",
             (const void *)numbering, cs_numbering_type_name[numbering->type],
             numbering->vector_size,
             numbering->n_threads, numbering->n_groups,
             numbering->n_no_adj_halo_groups,
             (long)(numbering->n_no_adj_halo_elts));

  if (numbering->group_index != NULL) {

    bft_printf("\n  group start index:\n"
               "\n    group_id thread_id (id) start_index\n");

    for (i = 0; i < numbering->n_groups; i++) {
      for (j = 0; j < numbering->n_threads; j++) {
        int k = i*numbering->n_threads + j;
        bft_printf("      %2d       %2d      %3d   %d\n",
                   i, j, k, (int)(numbering->group_index[k]));
      }
      bft_printf("                       %3d\n",
                 numbering->n_groups*numbering->n_threads);

    }
  }

  bft_printf("\n\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
