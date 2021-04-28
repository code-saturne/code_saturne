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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
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

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  double t_measure = 1.0;
  double test_sum = 0.0;

  /* Initialization and environment */

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

  /* Finalize */

  exit (EXIT_SUCCESS);
}
