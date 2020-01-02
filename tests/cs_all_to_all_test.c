/*============================================================================
 * Unit test for cs_all_to_all.c;
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

#include "cs_defs.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_all_to_all.h"

/*---------------------------------------------------------------------------*/

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
    sprintf (filename, "cs_all_to_all_test_out.%d", rank);
    f = fopen(filename, "w");
    assert(f != NULL);
  }

  return vfprintf(f, format, arg_ptr);
}

static int
_bft_printf_flush_proxy(void)
{
  return fflush(NULL);
}

/*----------------------------------------------------------------------------
 * Stop the code in case of error
 *----------------------------------------------------------------------------*/

static void
_bft_error_handler(const char  *filename,
                   int          line_num,
                   int          sys_err_code,
                   const char  *format,
                   va_list      arg_ptr)
{
  CS_UNUSED(filename);
  CS_UNUSED(line_num);

  bft_printf_flush();

  if (sys_err_code != 0)
    fprintf(stderr, "\nSystem error: %s\n", strerror(sys_err_code));

  vfprintf(stderr, format, arg_ptr);
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  char mem_trace_name[32];
  int size = 1;
  int rank = 0;

  int *dest_rank = NULL;

#if defined(HAVE_MPI)

  /* Initialization */

  cs_base_mpi_init(&argc, &argv);

  if (cs_glob_mpi_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(cs_glob_mpi_comm, &rank);
    MPI_Comm_size(cs_glob_mpi_comm, &size);
  }

  if (size < 1)
    return 0;

#endif /* (HAVE_MPI) */

  bft_error_handler_set(_bft_error_handler);
  bft_printf_proxy_set(_bft_printf_proxy);
  bft_printf_flush_proxy_set(_bft_printf_flush_proxy);

  sprintf(mem_trace_name, "cs_all_to_all_test_mem.%d", rank);
  bft_mem_init(mem_trace_name);

  cs_all_to_all_type_t a2at[3] = {CS_ALL_TO_ALL_MPI_DEFAULT,
                                  CS_ALL_TO_ALL_CRYSTAL_ROUTER,
                                  CS_ALL_TO_ALL_CRYSTAL_ROUTER};

  int a2a_flags[3] = {0, 0, CS_ALL_TO_ALL_ORDER_BY_SRC_RANK};

  for (int test_id = 0; test_id < 3; test_id++) {

    cs_all_to_all_set_type(a2at[test_id]);

    int flags = a2a_flags[test_id];

    bft_printf("\n"
               "Using all-to-all type %d (flags %d)\n"
               "---------------------\n\n",
               (int)a2at[test_id], flags);

    /* Build test array */

    cs_lnum_t n_elts = 3 + rank%3;

    BFT_MALLOC(dest_rank, n_elts, int);

    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      int _rank = rank + ii%5 - 2;
      if (_rank < 0) _rank = 0;
      if (_rank >= size) _rank = size-1;
      dest_rank[ii] = _rank;
    }

    cs_lnum_t *src_index = NULL;
    BFT_MALLOC(src_index, n_elts + 1, cs_lnum_t);

    src_index[0] = 0;
    for (cs_lnum_t ii = 0; ii < n_elts; ii++)
      src_index[ii+1] = src_index[ii] + 2 + ii%2;

    cs_gnum_t *src_val = NULL;
    BFT_MALLOC(src_val, src_index[n_elts], cs_gnum_t);
    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      bft_printf("%d -> rank %d :", ii, dest_rank[ii]);
      src_val[src_index[ii]] = ii;
      src_val[src_index[ii]+1] = rank;
      for (cs_lnum_t jj = src_index[ii] + 2;
           jj < src_index[ii+1];
           jj++)
        src_val[jj] = jj;
      for (cs_lnum_t jj = src_index[ii];
           jj < src_index[ii+1];
           jj++)
        bft_printf(" %llu", (unsigned long long)src_val[jj]);
      bft_printf("\n");
    }

    cs_all_to_all_t *d = cs_all_to_all_create(n_elts,
                                              flags,
                                              NULL, /* dest_id */
                                              dest_rank,
                                              cs_glob_mpi_comm);


    cs_lnum_t *dest_index = cs_all_to_all_copy_index(d,
                                                     false, /* reverse */
                                                     src_index,
                                                     NULL);

    cs_gnum_t *dest_val = cs_all_to_all_copy_indexed(d,
                                                     CS_GNUM_TYPE,
                                                     false, /* reverse */
                                                     src_index,
                                                     src_val,
                                                     dest_index,
                                                     NULL);

    cs_lnum_t n_elts_dest = cs_all_to_all_n_elts_dest(d);

    for (cs_lnum_t ii = 0; ii < n_elts_dest; ii++) {
      bft_printf("r %d -> (%d - %d) :", ii, dest_index[ii], dest_index[ii+1]);
      for (cs_lnum_t jj = dest_index[ii];
           jj < dest_index[ii+1];
           jj++)
        bft_printf(" %llu", (unsigned long long)dest_val[jj]);
      bft_printf("\n");
    }

    bft_printf("\nPrepare reverse\n\n");

    cs_gnum_t *reverse_val = NULL;
    BFT_MALLOC(reverse_val, dest_index[n_elts_dest] + n_elts_dest, cs_gnum_t);

    /* insert one value per element for return */
    cs_lnum_t s_id = 0;
    for (cs_lnum_t ii = 0; ii < n_elts_dest; ii++) {
      cs_lnum_t n_sub = dest_index[ii+1] - s_id;
      for (cs_lnum_t jj = 0; jj < n_sub; jj++)
        reverse_val[dest_index[ii] + jj] = dest_val[s_id + jj];
      reverse_val[dest_index[ii] + n_sub] = ii+1;
      s_id = dest_index[ii+1];
      dest_index[ii+1] += ii+1;
      bft_printf("%d -> (%d - %d) :", ii, dest_index[ii], dest_index[ii+1]);
      for (cs_lnum_t jj = dest_index[ii];
           jj < dest_index[ii+1];
           jj++)
        bft_printf(" %llu", (unsigned long long)reverse_val[jj]);
      bft_printf("\n");
    }

    BFT_FREE(dest_val);

    cs_all_to_all_copy_index(d,
                             true, /* reverse */
                             dest_index,
                             src_index);

    cs_gnum_t *ret_val = cs_all_to_all_copy_indexed(d,
                                                    CS_GNUM_TYPE,
                                                    true, /* reverse */
                                                    dest_index,
                                                    reverse_val,
                                                    src_index,
                                                    NULL);

    cs_all_to_all_destroy(&d);

    BFT_FREE(reverse_val);

    for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
      bft_printf("r %d -> (%d - %d) :", ii, src_index[ii], src_index[ii+1]);
      for (cs_lnum_t jj = src_index[ii];
           jj < src_index[ii+1];
           jj++)
        bft_printf(" %llu", (unsigned long long)ret_val[jj]);
      bft_printf("\n");
    }

    BFT_FREE(dest_index);
    BFT_FREE(ret_val);

    BFT_FREE(dest_rank);

    BFT_FREE(src_index);
    BFT_FREE(src_val);

  }

  bft_mem_end();

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Finalize();
  }
#endif

  exit (EXIT_SUCCESS);
}
