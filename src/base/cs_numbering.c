/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Numbering information for vectorization or multithreading
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

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

const char  *cs_numbering_type_name[] = {N_("vectorization"),
                                         N_("threads")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a numbering information structure in case of vectorization.
 *
 * parameters:
 *   vector_size <-- vector size used for this vectorization
 *
 * returns:
 *   pointer to created cs_numbering_t structure
 *---------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_vectorized(int  vector_size)
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_VECTORIZE;

  numbering->vector_size = vector_size;

  numbering->n_threads = 1;
  numbering->n_groups = 1;

  numbering->group_index = NULL;

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
cs_numbering_create_threaded(int         n_threads,
                             int         n_groups,
                             fvm_lnum_t  group_index[])
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_THREADS;

  numbering->vector_size = 1;

  numbering->n_threads = n_threads;
  numbering->n_groups = n_groups;

  BFT_MALLOC(numbering->group_index, n_threads*2*n_groups, fvm_lnum_t);

  memcpy(numbering->group_index,
         group_index,
         (n_threads*2*n_groups) * sizeof(fvm_lnum_t));

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

    if (_n->group_index != NULL)
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

  bft_printf("\n  Numbering:         %p\n"
             "  type:           %s\n"
             "  vector_size:    %d\n"
             "  n_threads:      %d\n"
             "  n_groups:       %d\n",
             numbering, cs_numbering_type_name[numbering->type],
             numbering->n_threads, numbering->n_groups);

  if (numbering->group_index != NULL) {

    bft_printf("\n  group start index:\n"
               "\n    group_id thread_id (id) start_index\n");

    for (i = 0; i < numbering->n_groups; i++) {
      for (j = 0; j < numbering->n_threads; j++) {
        int k = i*numbering->n_groups + j;
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
