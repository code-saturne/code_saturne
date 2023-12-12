/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

#include "cs_system_info.h"

#include "cs_defs.h"
#include "cs_math.h"
#include "cs_timer.h"

#include "cs_base_accel.h"

extern "C" {
  void
  main_cuda(void);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * False print of a message to standard output for discarded logs
 *----------------------------------------------------------------------------*/

static int
_bft_printf_null(const char  *format,
                 va_list      arg_ptr)
{
  CS_UNUSED(format);
  CS_UNUSED(arg_ptr);

  return 0;
}

/*----------------------------------------------------------------------------
 * Analysis of environment variables to determine
 * if we require MPI, and initialization if necessary.
 *----------------------------------------------------------------------------*/

static void
_mpi_init(void)
{
  int flag = 0;

  MPI_Initialized(&flag);

  if (!flag) {
#if defined(MPI_VERSION) && (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
    int mpi_threads;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &mpi_threads);
#else
    MPI_Init(NULL, NULL);
#endif
  }

  cs_glob_mpi_comm = MPI_COMM_WORLD;
  MPI_Comm_size(cs_glob_mpi_comm, &cs_glob_n_ranks);
  MPI_Comm_rank(cs_glob_mpi_comm, &cs_glob_rank_id);

  if (cs_glob_rank_id > 0)
    bft_printf_proxy_set(_bft_printf_null);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * OpenMP offload test.
 *----------------------------------------------------------------------------*/

static void
_omp_target_test(void)
{
#if defined(HAVE_OPENMP_TARGET)
  int m = 10, n = 500;
  double a[n][m], b[n][m], c[n][m];

  int n_devices = omp_get_num_devices();

  printf("Number of OpenMP target devices: %d\n"
         "Selected OpenMP target device:   %d\n",
         n_devices, cs_get_device_id());

  #pragma omp target
  {
    if (omp_is_initial_device())
      printf("  Running on host\n");
    else
      printf("  Running on device\n");
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
       a[i][j] = 1;
       b[i][j] = 2;
    }
  }

  #pragma omp target data map(to: a, b) map(from: c)
  {
    #pragma omp target teams distribute
    for (int i = 0; i < n; i++) {
      #pragma omp parallel for
      for (int j = 0; j < m; j++)
        c[i][j] = a[i][j] + b[i][j];
    }
  }

  printf("OpenMP target test result: %g\n", c[n-1][m-1]);
#endif
}

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  /* Initialization and environment */

#if defined(HAVE_MPI)
  _mpi_init();
#endif

  if (getenv("CS_MEM_LOG") != NULL) {
    char mem_log_file_name[128];
    int r_id = CS_MAX(cs_glob_rank_id, 0);
    snprintf(mem_log_file_name, 127, "%s.%d",
             getenv("CS_MEM_LOG"), r_id);
    bft_mem_init(mem_log_file_name);
  }
  else
    bft_mem_init(NULL);

  (void)cs_timer_wtime();

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif

  bft_printf("\n");

  /* Allocation tests */
  /*------------------*/

#if defined(HAVE_OPENMP_TARGET)

  cs_omp_target_select_default_device();  /* Initialize device id */

#endif

  cs_real_t *a0, *a1, *a2;
  CS_MALLOC_HD(a0, 100, cs_real_t, CS_ALLOC_HOST);

  if (cs_get_device_id() > -1) {
    CS_MALLOC_HD(a1, 100, cs_real_t, CS_ALLOC_HOST_DEVICE);
    CS_MALLOC_HD(a2, 100, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);
  }
  else {
    CS_MALLOC_HD(a1, 100, cs_real_t, CS_ALLOC_HOST);
    CS_MALLOC_HD(a2, 100, cs_real_t, CS_ALLOC_HOST);
  }

  bft_printf("Number of current allocations: %d\n", cs_get_n_allocations_hd());

  CS_FREE_HD(a0);
  CS_FREE_HD(a1);
  CS_FREE_HD(a2);

  bft_printf("Number of current allocations: %d\n", cs_get_n_allocations_hd());

  /* OpenMP tests */

  _omp_target_test();

  /* CUDA tests */

#if defined(HAVE_CUDA)
  main_cuda();
#endif

  /* Finalize */

  bft_mem_end();

  /* Finalize */

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Finalize();
  }
#endif /* HAVE_MPI */

  exit (EXIT_SUCCESS);
}
