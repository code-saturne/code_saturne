/*============================================================================
 * Unit test for fvm_interface.c;
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2007-2009  EDF

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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "fvm_interface.h"
#include "fvm_parall.h"
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

static fvm_interface_set_t *
_periodic_is(int                       comm_size,
             int                       comm_rank,
             fvm_gnum_t                n_side,
             int                       n_periodic_lists,
             const fvm_periodicity_t  *perio)
{
  fvm_gnum_t ii, jj, kk, couple_id;

  fvm_gnum_t n_max = n_side * n_side * n_side;
  fvm_lnum_t n_elements = ((n_max+1)*4) / (3*comm_size);
  int periodicity_num[3] = {1, 2, 3};

  fvm_lnum_t *n_periodic_couples = NULL;
  fvm_gnum_t **periodic_couples = NULL;
  fvm_gnum_t *global_number = NULL;

  fvm_interface_set_t *ifset = NULL;

  if (comm_size > 1) {

    BFT_MALLOC(global_number, n_elements, fvm_gnum_t);

    for (ii = 0; ii < (fvm_gnum_t)n_elements; ii++) {
      global_number[ii] = (n_elements * 3 / 4) * comm_rank + ii + 1;
    }
    if (comm_rank == comm_size -1 && global_number[n_elements - 1] < n_max)
      global_number[n_elements - 1] = n_max;

  }

  BFT_MALLOC(n_periodic_couples, n_periodic_lists, fvm_lnum_t);
  BFT_MALLOC(periodic_couples, n_periodic_lists, fvm_gnum_t *);

  for (ii = 0; ii < (fvm_gnum_t)n_periodic_lists; ii++) {

    n_periodic_couples[ii] = 0;
    periodic_couples[ii] = NULL;

    if (comm_rank == 0 || comm_rank == 1) {

      n_periodic_couples[ii] = n_side*n_side;
      BFT_MALLOC(periodic_couples[ii], n_periodic_couples[ii]*2, fvm_gnum_t);

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
      BFT_REALLOC(periodic_couples[ii], n_periodic_couples[ii]*2, fvm_gnum_t);

    }

  }

  ifset = fvm_interface_set_create(n_elements,
                                   NULL,
                                   global_number,
                                   perio,
                                   n_periodic_lists,
                                   periodicity_num,
                                   n_periodic_couples,
                                   (const fvm_gnum_t **const)periodic_couples);

  for (ii = 0; ii < (fvm_gnum_t)n_periodic_lists; ii++)
    BFT_FREE(periodic_couples[ii]);
  BFT_FREE(periodic_couples);
  BFT_FREE(n_periodic_couples);

  if (comm_size > 1)
    BFT_FREE(global_number);

  return ifset;
}

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message sur la sortie standard
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    sprintf (filename, "fvm_interface_test_out.%d", rank);
    f = fopen(filename, "w");
    assert(f != NULL);
  }

  return vfprintf(f, format, arg_ptr);
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  char mem_trace_name[32];
  int size = 1;
  int rank = 0;
  fvm_interface_set_t *ifset = NULL;
  fvm_periodicity_t *perio = NULL;

#if defined(HAVE_MPI)

  int ii;

  MPI_Status status;

  int sync = 1;

  fvm_lnum_t n_elements = 20;
  fvm_gnum_t *global_number = NULL;

  /* Initialization */

  MPI_Init(&argc, &argv);

  fvm_parall_set_mpi_comm(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size > 1)
    sprintf(mem_trace_name, "fvm_interface_test_mem.%d", rank);
  else
    strcpy(mem_trace_name, "fvm_interface_test_mem");
  bft_mem_init(mem_trace_name);

  bft_printf_proxy_set(_bft_printf_proxy);

  /* Build arbitrary interface */

  BFT_MALLOC(global_number, n_elements, fvm_gnum_t);

  for (ii = 0; ii < n_elements; ii++) {
    global_number[ii] = (n_elements * 3 / 4) * rank + ii + 1;
  }

  ifset = fvm_interface_set_create(n_elements,
                                   NULL,
                                   global_number,
                                   NULL,
                                   0,
                                   NULL,
                                   NULL,
                                   NULL);

  BFT_FREE(global_number);

  /* Serialize dump of interfaces */

  if (rank > 0)
    MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);

  bft_printf("Interface on rank %d:\n\n", rank);

  fvm_interface_set_dump(ifset);

  if (rank < size - 1)
    MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);

  /* We are finished for this interface */

  ifset = fvm_interface_set_destroy(ifset);

#else

  strcpy(mem_trace_name, "fvm_interface_test_mem");
  bft_mem_init(mem_trace_name);

  bft_printf_proxy_set(_bft_printf_proxy);

#endif /* (HAVE_MPI) */

  /* Now build interface with periodicity */

  perio = _init_periodicity();

  /* Dump periodicity info for first rank (identical on all ranks) */

  if (rank == 0)
    fvm_periodicity_dump(perio);

  ifset = _periodic_is(size,
                       rank,
                       4, /* n vertices/cube side */
                       3, /* n_periodic_lists */
                       perio);

#if defined(HAVE_MPI)

  /* Serialize dump of interfaces */

  if (rank > 0)
    MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);

#endif /* (HAVE_MPI) */

  bft_printf("Periodic Interface on rank %d:\n\n", rank);

  fvm_interface_set_dump(ifset);

#if defined(HAVE_MPI)

  if (rank < size - 1)
    MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);

#endif /* (HAVE_MPI) */

  /* We are finished */

  ifset = fvm_interface_set_destroy(ifset);
  fvm_periodicity_destroy(perio);

  bft_mem_end();

#if defined(HAVE_MPI)

  MPI_Finalize();

#endif

  exit (EXIT_SUCCESS);
}
