/*============================================================================
 * Unit test for fvm_file.c;
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2007-2008  EDF

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_file.h>
#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "fvm_file.h"
#include "fvm_parall.h"

/*---------------------------------------------------------------------------*/

static void
_create_test_data(void)
{
  int i;
  bft_file_t *f;

  char header[80];
  char footer[80];
  int iarray[30];
  double farray[30];

  sprintf(header, "fvm test file");
  for (i = strlen(header); i < 80; i++)
    header[i] = '\0';

  sprintf(footer, "fvm test file end");
  for (i = strlen(footer); i < 80; i++)
    footer[i] = '\0';

  for (i = 0; i < 30; i++)
    iarray[i] = i+1;
  for (i = 0; i < 30; i++)
    farray[i] = i+1;

  f= bft_file_open("file_test_data",
                   BFT_FILE_MODE_WRITE,
                   BFT_FILE_TYPE_BINARY);

  bft_file_set_big_endian(f);

  bft_file_write(header, 1, 80, f);
  bft_file_write(iarray, sizeof(int), 30, f);
  bft_file_write(farray, sizeof(double), 30, f);

  bft_file_write(footer, 1, 80, f);

  f = bft_file_free(f);
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  char mem_trace_name[32];
  char output_file_name[32];
  char buf[80];
  int ibuf[30];
  double dbuf[30];
  fvm_gnum_t i;
  int test_id;
  int size = 1;
  int rank = 0;
  size_t retval = 0;
  fvm_gnum_t block_start, block_end;
  fvm_gnum_t block_start_2, block_end_2;
  fvm_file_off_t off1 = -1, off2 = -1;
  fvm_file_t *f = NULL;


#if defined(FVM_HAVE_MPI_IO)
  const int n_hints = 3;
  const int hints[3] = {FVM_FILE_NO_MPI_IO,
                        FVM_FILE_EXPLICIT_OFFSETS,
                        FVM_FILE_INDIVIDUAL_POINTERS};
#else
  const int n_hints = 1;
  const int hints[1] = {FVM_FILE_NO_MPI_IO};
#endif

#if defined(FVM_HAVE_MPI)

  MPI_Status status;
  int sync = 1;

  /* Initialization */

  MPI_Init(&argc, &argv);

  fvm_parall_set_mpi_comm(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  block_start = rank * (30./size) + 1;
  block_end   = (rank + 1) * (30./size) + 1;
  if (rank == size - 1)
    block_end = 31;

  block_start_2 = rank * (15./size) + 1;
  block_end_2   = (rank + 1) * (15./size) + 1;
  if (rank == size - 1)
    block_end_2 = 16;

#else

  block_start = 1;
  block_end = 31;

  block_start_2 = 1;
  block_end_2 = 16;

#endif /* (FVM_HAVE_MPI) */

  if (size > 1)
    sprintf(mem_trace_name, "fvm_file_test_mem.%d", rank);
  else
    strcpy(mem_trace_name, "fvm_file_test_mem");
  bft_mem_init(mem_trace_name);

  if (rank == 0)
    _create_test_data();

  /* Loop on tests */

  for (test_id = 0; test_id < n_hints; test_id++) {

    if (rank == 0)
      bft_printf("Running test: %d\n"
                 "-------------\n\n", test_id);

    sprintf(output_file_name, "output_data_%d", test_id+1);

    /* Read and seek/set tests */
    /*-------------------------*/

#if defined(FVM_HAVE_MPI)

    f = fvm_file_open("file_test_data",
                      FVM_FILE_MODE_READ,
                      hints[test_id],
                      MPI_COMM_WORLD);

#else

    f = fvm_file_open("file_test_data",
                      FVM_FILE_MODE_READ,
                      hints[test_id]);

#endif /* (FVM_HAVE_MPI) */

    fvm_file_set_big_endian(f);

    fvm_file_dump(f);

    retval = fvm_file_read_global(f, buf, 1, 80);

    bft_printf("rank %d, readbuf = %s (returned %d)\n\n",
               rank, buf, (int)retval);

    for (i = 0; i < 30; i++)
      ibuf[i] = 0;

    retval = fvm_file_read_block(f, ibuf, sizeof(int), 1,
                                 block_start, block_end);

#if defined(FVM_HAVE_MPI) /* Serialize dump */
    if (rank > 0)
      MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
#endif

    bft_printf("\nRead by rank %d (returned %d):\n\n", rank, (int)retval);
    for (i = block_start; i < block_end; i++)
      bft_printf("  ival[%d] = %d\n", (int)i, (int)(ibuf[i-block_start]));

#if defined(FVM_HAVE_MPI)
    if (rank < size - 1)
      MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
#endif

    off1 = fvm_file_tell(f);

    retval = fvm_file_read_block(f, dbuf, sizeof(double), 2,
                                 block_start_2, block_end_2);

    off2 = fvm_file_tell(f);

#if defined(FVM_HAVE_MPI) /* Serialize dump */
    if (rank > 0)
      MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
#endif

    bft_printf("\nOffsets saved by rank %d: %ld, %ld\n\n",
               rank, (long)off1, (long)off2);

    bft_printf("\nRead by rank %d (returned %d):\n\n", rank, (int)retval);
    for (i = block_start_2; i < block_end_2; i++) {
      bft_printf("  dval[%d] = %f\n", (int)(i*2 - 1),
                 (dbuf[(i-block_start_2)*2]));
      bft_printf("  dval[%d] = %f\n", (int)(i*2),
                 (dbuf[(i-block_start_2)*2+1]));
    }

#if defined(FVM_HAVE_MPI)
    if (rank < size - 1)
      MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
#endif

#if defined(FVM_HAVE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    bft_printf("barrier passed by rank %d\n", rank);

    retval = fvm_file_read_global(f, buf, 1, 80);
    bft_printf("rank %d, buf = %s (returned %d)\n", rank, buf, (int)retval);

    /* Test seek by re-reading at saved offset */

    fvm_file_seek(f, off1, FVM_FILE_SEEK_SET);

    memset(dbuf, 0, (block_end - block_start)*sizeof(double));
    retval = fvm_file_read_block(f, dbuf, sizeof(double), 1,
                                 block_start, block_end);

#if defined(FVM_HAVE_MPI) /* Serialize dump */
    if (rank > 0)
      MPI_Recv(&sync, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
#endif

    bft_printf("\nRe-read by rank %d (returned %d):\n\n", rank, (int)retval);
    for (i = block_start; i < block_end; i++)
      bft_printf("  dval[%d] = %f\n", (int)i, (dbuf[i-block_start]));

#if defined(FVM_HAVE_MPI)
    if (rank < size - 1)
      MPI_Send(&sync, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
#endif

    fvm_file_seek(f, off2, FVM_FILE_SEEK_SET);

    retval = fvm_file_read_global(f, buf, 1, 80);
    bft_printf("rank %d, re-read buf = %s (returned %d)\n",
               rank, buf, (int)retval);

    f = fvm_file_free(f);

    /* Write tests */
    /*-------------*/

#if defined(FVM_HAVE_MPI)

    f = fvm_file_open(output_file_name,
                      FVM_FILE_MODE_WRITE,
                      hints[test_id],
                      MPI_COMM_WORLD);

#else

    f = fvm_file_open(output_file_name,
                      FVM_FILE_MODE_WRITE,
                      hints[test_id]);

#endif /* (FVM_HAVE_MPI) */

    fvm_file_set_big_endian(f);

    fvm_file_dump(f);

    sprintf(buf, "fvm test file");
    for (i = strlen(buf); i < 80; i++)
      buf[i] = '\0';

    retval = fvm_file_write_global(f, buf, 1, 80);

    bft_printf("rank %d, wrote %d global values.\n", rank, (int)retval);

    for (i = block_start_2; i < block_end_2; i++) {
      ibuf[(i-block_start_2)*2] = i*2 - 1;
      ibuf[(i-block_start_2)*2 + 1] = i*2;
    }
    for (i = block_start; i < block_end; i++)
      dbuf[i-block_start] = i;

    retval = fvm_file_write_block(f, ibuf, sizeof(int), 2,
                                  block_start_2, block_end_2);

    bft_printf("rank %d, wrote %d block values.\n", rank, (int)retval);

    retval = fvm_file_write_block_buffer(f, dbuf, sizeof(double), 1,
                                         block_start, block_end);

    bft_printf("rank %d, wrote %d block (buffer) values.\n",
               rank, (int)retval);

    sprintf(buf, "fvm test file end");
    for (i = strlen(buf); i < 80; i++)
      buf[i] = '\0';

    retval = fvm_file_write_global(f, buf, 1, 80);

    bft_printf("rank %d, wrote %d global values.\n", rank, (int)retval);

    f = fvm_file_free(f);
  }

  /* We are finished */

  bft_mem_end();

#if defined(FVM_HAVE_MPI)

  MPI_Finalize();

#endif

  exit (EXIT_SUCCESS);
}
