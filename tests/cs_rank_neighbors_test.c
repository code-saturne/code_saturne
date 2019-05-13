/*============================================================================
 * Unit test for cs_rank_neighbors.c;
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_rank_neighbors.h"

/*---------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Print message on standard output
 *----------------------------------------------------------------------------*/

static int _bft_printf_proxy
(
 const char     *const format,
       va_list         arg_ptr
)
{
  static FILE *f = NULL;

  if (f == NULL) {
    char filename[64];
    int rank = 0;
#if defined(HAVE_MPI)
    if (cs_glob_mpi_comm != MPI_COMM_NULL)
      MPI_Comm_rank(cs_glob_mpi_comm, &rank);
#endif
    sprintf (filename, "cs_rank_neighbors_test_out.%d", rank);
    f = fopen(filename, "w");
    assert(f != NULL);
  }

  return vfprintf(f, format, arg_ptr);
}

/*----------------------------------------------------------------------------
 * Stop the code in case of error
 *----------------------------------------------------------------------------*/

static void
_bft_error_handler(const char  *filename,
                   int          line_no,
                   int          code_err_sys,
                   const char  *format,
                   va_list      arg_ptr)
{
  CS_UNUSED(filename);
  CS_UNUSED(line_no);

  bft_printf_flush();

  if (code_err_sys != 0)
    fprintf(stderr, "\nSystem error: %s\n", strerror(code_err_sys));

  vfprintf(stderr, format, arg_ptr);
}

/*----------------------------------------------------------------------------
 * Create a test rank neighbors structure
 *----------------------------------------------------------------------------*/

static cs_rank_neighbors_t *
_build_rank_neighbors(int rank,
                      int size,
                      int test_to_index)
{
  size_t n_elts = 10000;
  int *elt_rank;
  BFT_MALLOC(elt_rank, n_elts, int);

  int delta = sqrt(size);
  if (delta < 2)
    delta = 2;
  for (size_t i = 0; i < n_elts; i++)
    elt_rank[i] = (rank + (i/10)%delta) % size;

  cs_rank_neighbors_t *n = cs_rank_neighbors_create(n_elts,
                                                    elt_rank);

  if (test_to_index) {
    cs_rank_neighbors_to_index(n, n_elts, elt_rank, elt_rank);
    for (size_t i = 0; i < n_elts; i++) {
      assert(elt_rank[i] >= 0);
      assert(elt_rank[i] < n->size);
    }
    cs_lnum_t *elt_rank_count;
    BFT_MALLOC(elt_rank_count, n->size, cs_lnum_t);
    cs_rank_neighbors_count(n,
                            n_elts,
                            elt_rank,
                            elt_rank_count);
    bft_printf("Rank neighbors and counts on rank %d:\n", rank);
    for (int i = 0; i < n->size; i++)
      bft_printf("  %d: %d %d\n", i, n->rank[i], (int)(elt_rank_count[i]));
    BFT_FREE(elt_rank_count);
    bft_printf("\n");
  }

  BFT_FREE(elt_rank);

  return n;
}

#endif /* HAVE_MPI */

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  char mem_trace_name[32];
  int size = 1;
  int rank = 0;

#if defined(HAVE_MPI)

  /* Initialization */

  cs_base_mpi_init(&argc, &argv);

  if (cs_glob_mpi_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(cs_glob_mpi_comm, &rank);
    MPI_Comm_size(cs_glob_mpi_comm, &size);
  }

  bft_error_handler_set(_bft_error_handler);

  if (size > 1)
    sprintf(mem_trace_name, "cs_rank_neighbors_test_mem.%d", rank);
  else
    strcpy(mem_trace_name, "cs_rank_neighbors_test_mem");
  bft_mem_init(mem_trace_name);
  bft_printf_proxy_set(_bft_printf_proxy);

  cs_rank_neighbors_t *n = _build_rank_neighbors(rank, size, 1);

  bft_printf("Rank neighbors on rank %d:\n", rank);
  for (int i = 0; i < n->size; i++)
    bft_printf("  %d: %d\n", i, n->rank[i]);

  cs_rank_neighbors_symmetrize(n, cs_glob_mpi_comm);

  bft_printf("\nSymmetrized rank neighbors on rank %d:\n", rank);
  for (int i = 0; i < n->size; i++)
    bft_printf("%d: %d\n", i, n->rank[i]);

  /* Other methods */

  for (cs_rank_neighbors_exchange_t t = CS_RANK_NEIGHBORS_NBX;
       t <= CS_RANK_NEIGHBORS_CRYSTAL_ROUTER;
       t++) {

    cs_rank_neighbors_set_exchange_type(t);
    if (cs_rank_neighbors_get_exchange_type() != t) {
      bft_printf("%s not available\n",
                 cs_rank_neighbors_exchange_name[t]);
      continue;
    }

    cs_rank_neighbors_t *nc = _build_rank_neighbors(rank, size, 0);

    cs_rank_neighbors_symmetrize(nc, cs_glob_mpi_comm);

    assert(n->size == nc->size);
    for (int i = 0; i < n->size; i++) {
      assert(n->rank[i] == nc->rank[i]);
    }

    cs_rank_neighbors_destroy(&nc);

  }

  cs_rank_neighbors_destroy(&n);

  cs_rank_neighbors_set_exchange_type(CS_RANK_NEIGHBORS_PEX);

  n = _build_rank_neighbors(rank, size, 0);

  cs_rank_neighbors_t  *nr;
  cs_lnum_t            *send_count, *recv_count;

  BFT_MALLOC(send_count, n->size, cs_lnum_t);

  for (int i = 0; i < n->size; i++)
    send_count[i] = 5 + i%4;

  cs_rank_neighbors_sync_count(n, &nr, send_count, &recv_count,
                               cs_glob_mpi_comm);

  bft_printf("\nReceived on rank %d from:\n", rank);
  for (int i = 0; i < nr->size; i++)
    bft_printf("  %d: %d (%d)\n", i, nr->rank[i], (int)(recv_count[i]));

  /* Other methods */

  for (cs_rank_neighbors_exchange_t t = CS_RANK_NEIGHBORS_NBX;
       t <= CS_RANK_NEIGHBORS_CRYSTAL_ROUTER;
       t++) {

    cs_rank_neighbors_t  *nr1;
    cs_lnum_t            *recv_count1;

    cs_rank_neighbors_set_exchange_type(t);
    if (cs_rank_neighbors_get_exchange_type() != t) {
      bft_printf("%s not available\n",
                 cs_rank_neighbors_exchange_name[t]);
      continue;
    }

    cs_rank_neighbors_sync_count(n, &nr1, send_count, &recv_count1,
                                 cs_glob_mpi_comm);

    assert(nr->size == nr1->size);
    for (int i = 0; i < nr->size; i++) {
      assert(nr->rank[i] == nr1->rank[i]);
      assert(recv_count[i] == recv_count1[i]);
    }

    BFT_FREE(recv_count1);
    cs_rank_neighbors_destroy(&nr1);

  }

  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  cs_rank_neighbors_destroy(&nr);
  cs_rank_neighbors_destroy(&n);

  bft_mem_end();

  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Finalize();
  }

#else

  bft_printf("rank neighborhoods only make sense for MPI\n");

#endif

  exit (EXIT_SUCCESS);
}
