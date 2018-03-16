/*============================================================================
 * Unit test for ple_coupling.c;
 *============================================================================*/

/*
  This file is part of the "Distributed Projection and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2018  EDF S.A.

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

  int rank, i, j;

  int n_apps = 1;
  int app_id = -1;
  int sync_flag = PLE_COUPLING_INIT;
  MPI_Comm app_comm = MPI_COMM_NULL;
  MPI_Comm intracomm = MPI_COMM_NULL;
  ple_coupling_mpi_set_t *s = NULL;
  ple_coupling_mpi_set_info_t info[4];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank < 2)
    app_id = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, "Code_A");
  else if (rank < 3)
    app_id = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, "Code_B:case_b");
  else if (rank < 6)
    app_id = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, "Code_C");
  else
    app_id = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, "Code_D;case_d");

  if (app_id < 0)
    app_id = 0;

  MPI_Comm_split(MPI_COMM_WORLD, app_id, rank, &app_comm);

  switch(app_id) {
  case 0:
    sync_flag = PLE_COUPLING_NO_SYNC;
    s = ple_coupling_mpi_set_create(sync_flag, "Code_A", NULL,
                                    MPI_COMM_WORLD, app_comm);
    break;
  case 1:
    s = ple_coupling_mpi_set_create(sync_flag, "Code_B", "case b",
                                    MPI_COMM_WORLD, app_comm);
    break;
  case 2:
    s = ple_coupling_mpi_set_create(sync_flag, "Code_C", NULL,
                                    MPI_COMM_WORLD, app_comm);
    break;
  default:
    s = ple_coupling_mpi_set_create(sync_flag, "Code_D", "case d",
                                    MPI_COMM_WORLD, app_comm);
    break;
  }

  ple_coupling_mpi_set_dump(s);

  for (i = 0; i < 4; i++) {

    const int *status = NULL;
    const double *ts = NULL;

    switch(app_id) {
    case 0:
      break;
    case 1:
      ple_coupling_mpi_set_synchronize(s, PLE_COUPLING_NEW_ITERATION, 0.1);
      break;
    case 2:
      if (i == 1)
        ple_coupling_mpi_set_synchronize(s,
                                         (  PLE_COUPLING_NEW_ITERATION
                                          | PLE_COUPLING_LAST),
                                         0.1);
      else
        ple_coupling_mpi_set_synchronize(s,
                                         PLE_COUPLING_NEW_ITERATION,
                                         0.2);
      break;
    default:
      ple_coupling_mpi_set_synchronize(s, PLE_COUPLING_NEW_ITERATION, 0.3);
      break;
    }

    if (app_id != 0) {
      status = ple_coupling_mpi_set_get_status(s);
      ts = ple_coupling_mpi_set_get_timestep(s);

      ple_printf("\n");
      for (j = 0; j < ple_coupling_mpi_set_n_apps(s); j++)
        ple_printf("app %d status %d, time step %f\n", j, status[j], ts[j]);
    }
  }

  n_apps = ple_coupling_mpi_set_n_apps(s);

  for (i = 0; i < n_apps; i++)
    info[i] = ple_coupling_mpi_set_get_info(s, i);

  ple_coupling_mpi_set_destroy(&s);

  if (n_apps > 1 && app_id < 2) {

    int local_range[2],  distant_range[2];
    int dist_app_id = (app_id + 1) % 2;
    int dist_root_rank = info[dist_app_id].root_rank;

    local_range[0] = info[app_id].root_rank;
    local_range[1] = local_range[0] + info[app_id].n_ranks;
    distant_range[0] = info[dist_app_id].root_rank;
    distant_range[1] = distant_range[0] + info[dist_app_id].n_ranks;

    ple_printf("\nBuilding intracomm with:\n"
               "  dist_root_rank = %d\n"
               "  local_range =   [%d, %d]\n"
               "  distant_range = [%d, %d]\n",
               dist_root_rank,
               local_range[0], local_range[1],
               distant_range[0], distant_range[1]);

    ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                      app_comm,
                                      dist_root_rank,
                                      &intracomm,
                                      local_range,
                                      distant_range);

  }

  MPI_Finalize();

#endif /* (PLE_HAVE_MPI) */

  exit (EXIT_SUCCESS);
}
