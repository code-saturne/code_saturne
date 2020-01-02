/*============================================================================
 * Unit test for cs_interface.c;
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
#include "cs_interface.h"
#include "fvm_periodicity.h"

/*---------------------------------------------------------------------------*/

static fvm_periodicity_t *
_init_periodicity(void)
{
  double t1[3] = {1., 0., 0.};
  double t2[3] = {0., 1., 0.};
  double t3[3] = {0., 0., 1.};

  /* double t4[3] = {1.02, 0., 0.}; */
  /* double a1[3] = {0., 0., 1.}; */
  /* double i1[3] = {0., 0., 0.}; */

  fvm_periodicity_t *p = NULL;

  p = fvm_periodicity_create(0.1);

  fvm_periodicity_add_translation(p, 1, t1);
  fvm_periodicity_add_translation(p, 2, t2);
  fvm_periodicity_add_translation(p, 3, t3);

  /* fvm_periodicity_add_rotation(p, 2, 90., a1, i1); */
  /* fvm_periodicity_add_translation(p, 3, t3); */

  fvm_periodicity_combine(p, 1);

  return p;
}

static void
_periodic_is_test(int                       ordered_gnum,
                  int                       comm_size,
                  int                       comm_rank,
                  cs_gnum_t                 n_side,
                  int                       n_periodic_lists,
                  const fvm_periodicity_t  *perio)
{
  cs_lnum_t ii, couple_id;
  cs_gnum_t jj, kk;

  cs_gnum_t n_3 = n_side * n_side * n_side;
  cs_lnum_t n_elements = ((n_3+1)*4) / (3*comm_size);
  cs_gnum_t n_max = (n_elements * 3 / 4) * (comm_size-1) + n_elements;
  int periodicity_num[3] = {1, 2, 3};

  cs_lnum_t *n_periodic_couples = NULL;
  cs_gnum_t **periodic_couples = NULL;
  cs_gnum_t *global_number = NULL;
  cs_real_t *var = NULL;

  cs_interface_set_t *ifset = NULL;

  if (comm_size > 1) {

    BFT_MALLOC(global_number, n_elements, cs_gnum_t);

    for (ii = 0; ii < n_elements; ii++)
      global_number[ii] = (n_elements * 3 / 4) * comm_rank + ii + 1;

    if (comm_rank == comm_size -1 && global_number[n_elements - 1] < n_max)
      global_number[n_elements - 1] = n_max;

    if (ordered_gnum == false) {
      for (ii = 0; ii < n_elements; ii++) {
        global_number[ii] = n_max + 1 - global_number[ii];
      }
    }

  }

  BFT_MALLOC(n_periodic_couples, n_periodic_lists, cs_lnum_t);
  BFT_MALLOC(periodic_couples, n_periodic_lists, cs_gnum_t *);

  for (ii = 0; ii < n_periodic_lists; ii++) {

    n_periodic_couples[ii] = 0;
    periodic_couples[ii] = NULL;

    if (comm_rank == 0 || comm_rank == 1) {

      n_periodic_couples[ii] = n_side*n_side;
      BFT_MALLOC(periodic_couples[ii], n_periodic_couples[ii]*2, cs_gnum_t);

      couple_id = 0;

      switch(ii) {
      case 0:
        for (jj = 0; jj < n_side*n_side; jj++) {
          periodic_couples[ii][couple_id*2] = (jj*n_side) + 1;
          periodic_couples[ii][couple_id*2+1]
            = periodic_couples[ii][couple_id*2] + (n_side-1);
          couple_id++;
        }
        break;
      case 1:
        for (jj = 0; jj < n_side; jj++) {
          for (kk = 0; kk < n_side; kk++) {
            periodic_couples[ii][couple_id*2] = jj*n_side*n_side + kk + 1;
            periodic_couples[ii][couple_id*2+1]
              = periodic_couples[ii][couple_id*2] + (n_side-1)*n_side;
            couple_id++;
          }
        }
        break;
      case 2:
        for (jj = 0; jj < n_side*n_side; jj++) {
          periodic_couples[ii][couple_id*2] = jj + 1;
          periodic_couples[ii][couple_id*2+1]
            = periodic_couples[ii][couple_id*2] + (n_side-1)*n_side*n_side;
          couple_id++;
        }
        break;
      default:
        break;
      }

      n_periodic_couples[ii] = couple_id;
      BFT_REALLOC(periodic_couples[ii], n_periodic_couples[ii]*2, cs_gnum_t);

    }

    if (ordered_gnum == false) {
      for (couple_id = 0; couple_id < n_periodic_couples[ii]; couple_id++) {
        periodic_couples[ii][couple_id*2]
          = n_max + 1 - periodic_couples[ii][couple_id*2];
        periodic_couples[ii][couple_id*2+1]
          = n_max + 1 - periodic_couples[ii][couple_id*2+1];
      }
    }

  }

  ifset = cs_interface_set_create(n_elements,
                                  NULL,
                                  global_number,
                                  perio,
                                  n_periodic_lists,
                                  periodicity_num,
                                  n_periodic_couples,
                                  (const cs_gnum_t **const)periodic_couples);

  for (ii = 0; ii < n_periodic_lists; ii++)
    BFT_FREE(periodic_couples[ii]);
  BFT_FREE(periodic_couples);
  BFT_FREE(n_periodic_couples);

  if (comm_size > 1)
    BFT_FREE(global_number);

  /* Now that interface set is created, dump info and run tests */

  bft_printf("Periodic Interface on rank %d:\n\n", comm_rank);

  cs_interface_set_add_match_ids(ifset);
  cs_interface_set_dump(ifset);
  cs_interface_set_free_match_ids(ifset);

  bft_printf("Interface set sum test:\n\n");

  BFT_MALLOC(var, n_elements*2, cs_real_t);

  for (ii = 0; ii < n_elements; ii++) {
    var[ii] = 1.0;
    var[ii + n_elements] = 2.0;
  }
  cs_interface_set_sum(ifset,
                       n_elements,
                       2,
                       false,
                       CS_REAL_TYPE,
                       var);
  for (ii = 0; ii < n_elements; ii++)
    bft_printf("%2d: %12.3f %12.3f\n", (int)ii, var[ii], var[ii + n_elements]);
  BFT_FREE(var);

  cs_interface_set_destroy(&ifset);
}

/*----------------------------------------------------------------------------
 * Print message on standard output
 *----------------------------------------------------------------------------*/

static int _bft_printf_proxy_o
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
    sprintf (filename, "cs_interface_test_o_out.%d", rank);
    f = fopen(filename, "w");
    assert(f != NULL);
  }

  return vfprintf(f, format, arg_ptr);
}

static int _bft_printf_proxy_u
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
    sprintf (filename, "cs_interface_test_u_out.%d", rank);
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
  int ii;
  char mem_trace_name[32];
  int size = 1;
  int rank = 0;
  cs_interface_set_t *ifset = NULL;
  fvm_periodicity_t *perio = NULL;

  cs_lnum_t n_elements = 20;

  cs_gnum_t *global_number = NULL;

#if defined(HAVE_MPI)

  /* Initialization */

  cs_base_mpi_init(&argc, &argv);

  if (cs_glob_mpi_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(cs_glob_mpi_comm, &rank);
    MPI_Comm_size(cs_glob_mpi_comm, &size);
  }

#endif /* (HAVE_MPI) */

  bft_error_handler_set(_bft_error_handler);

  if (size > 1)
    sprintf(mem_trace_name, "cs_interface_test_mem.%d", rank);
  else
    strcpy(mem_trace_name, "cs_interface_test_mem");
  bft_mem_init(mem_trace_name);

  for (int ordered = 1; ordered >= 0; ordered--) {

    if (ordered)
      bft_printf_proxy_set(_bft_printf_proxy_o);
    else
      bft_printf_proxy_set(_bft_printf_proxy_u);

    /* Build arbitrary interface */

    if (size > 1) {

      BFT_MALLOC(global_number, n_elements, cs_gnum_t);

      for (ii = 0; ii < n_elements; ii++)
        global_number[ii] = (n_elements * 3 / 4) * rank + ii + 1;

      if (ordered == 0) {
        cs_gnum_t gnum_max = (n_elements * 3 / 4) * (size - 1) + n_elements;
        for (ii = 0; ii < n_elements; ii++)
          global_number[ii] = gnum_max + 3 - global_number[ii];
      }

      ifset = cs_interface_set_create(n_elements,
                                      NULL,
                                      global_number,
                                      NULL,
                                      0,
                                      NULL,
                                      NULL,
                                      NULL);

      BFT_FREE(global_number);

      bft_printf("Interface on rank %d:\n\n", rank);

      cs_interface_set_add_match_ids(ifset);
      cs_interface_set_dump(ifset);
      cs_interface_set_free_match_ids(ifset);

      /* We are finished for this interface */

      cs_interface_set_destroy(&ifset);

    }

    /* Now build interface with periodicity */

    perio = _init_periodicity();

    /* Dump periodicity info for first rank (identical on all ranks) */

    if (rank == 0)
      fvm_periodicity_dump(perio);

    _periodic_is_test(ordered,
                      size,
                      rank,
                      4, /* n vertices/cube side */
                      3, /* n_periodic_lists */
                      perio);

    /* We are finished */

    fvm_periodicity_destroy(perio);

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
