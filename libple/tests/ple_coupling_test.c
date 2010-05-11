/*============================================================================
 * Unit test for ple_coupling.c;
 *============================================================================*/

/*
  This file is part of the "Distributed Projection and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2010  EDF

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

#include "ple_config_defs.h"
#include "ple_defs.h"

#if defined(PLE_HAVE_MPI)
#include <mpi.h>
#endif

#include "ple_coupling.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
#if defined(PLE_HAVE_MPI)

  int rank;

  int app_id = -1;
  MPI_Comm app_comm = MPI_COMM_NULL;
  MPI_Comm intracomm = MPI_COMM_NULL;
  ple_coupling_mpi_world_t *w = NULL;
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank < 2)
    app_id = 0;
  else if (rank < 3)
    app_id = 1;
  else if (rank < 6)
    app_id = 2;
  else
    app_id = 3;

  MPI_Comm_split(MPI_COMM_WORLD, app_id, rank, &app_comm);

  switch(app_id) {
  case 0:
    w = ple_coupling_mpi_world_create(app_id, "Code_A", NULL, app_comm);
    break;
  case 1:
    w = ple_coupling_mpi_world_create(app_id, "Code_B", "case b", app_comm);
    break;
  case 2:
    w = ple_coupling_mpi_world_create(app_id, "Code_C", NULL, app_comm);
    break;
  default:
    w = ple_coupling_mpi_world_create(app_id, "Code_D", "case d", app_comm);
    break;
  }

  ple_coupling_mpi_world_dump(w);

  ple_coupling_mpi_world_destroy(&w);

  if (app_id < 2) {
    int dist_root_rank = app_id == 0 ? 2 : 0;
    ple_coupling_mpi_intracomm_create(app_comm,
                                      dist_root_rank,
                                      &intracomm,
                                      local_range,
                                      distant_range);
  }

  MPI_Finalize();

#endif /* (PLE_HAVE_MPI) */

  exit (EXIT_SUCCESS);
}
