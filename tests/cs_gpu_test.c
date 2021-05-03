/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

/*----------------------------------------------------------------------------
 * OpenMP offload test.
 *----------------------------------------------------------------------------*/

static void
_omp_target_test(void)
{
  int m = 10, n = 500;
  double a[n][m], b[n][m], c[n][m];

#if defined(_OPENMP)

  int n_devices = omp_get_num_devices();
  printf("Number of OpenMP target devices %d\n", n_devices);

#endif

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
}

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  /* Initialization and environment */

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

  cs_real_t *a0, *a1, *a2;
  CS_MALLOC_HD(a0, 100, cs_real_t, CS_ALLOC_HOST);
  CS_MALLOC_HD(a1, 100, cs_real_t, CS_ALLOC_HOST_DEVICE);
  CS_MALLOC_HD(a2, 100, cs_real_t, CS_ALLOC_HOST_DEVICE_SHARED);

  bft_printf("Number of current allocations: %d\n", cs_get_n_allocations_hd());

  CS_FREE_HD(a0);
  CS_FREE_HD(a1);
  CS_FREE_HD(a2);

  bft_printf("Number of current allocations: %d\n", cs_get_n_allocations_hd());

  /* OpenMP tests */

  _omp_target_test();

  /* Finalize */

  bft_mem_end();

  exit (EXIT_SUCCESS);
}
