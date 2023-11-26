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

#include "cs_blas.h"
#include "cs_defs.h"
#include "cs_timer.h"

#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_assembler.h"
#include "cs_matrix_util.h"

#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/* Global matrix info */

cs_gnum_t  _n_g_vtx = 0;
cs_gnum_t  _vtx_range[2];

cs_lnum_t  _n_vtx = 0, _n_edges = 0;

cs_gnum_t  *_g_vtx_id = NULL;
cs_lnum_2_t  *_edges = NULL;

bool _full_connect = true;
bool _symmetric = true;

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Analysis of environment variables to determine
 * if we require MPI, and initialization if necessary.
 *----------------------------------------------------------------------------*/

static void
_mpi_init(void)
{
  int flag = 0;
  bool use_mpi = false;

#if defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  use_mpi = true;

#elif defined(MPICH2) || defined(MPICH)
  if (getenv("PMI_RANK") != NULL)
    use_mpi = true;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_COMM_WORLD_RANK") != NULL)    /* OpenMPI 1.3 and above */
    use_mpi = true;

#endif /* Tests for known MPI variants */

  /* If we have determined from known MPI environment variables
     of command line arguments that we are running under MPI,
     initialize MPI */

  if (use_mpi == true) {

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

  }

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * build base data based on rank id.
 *----------------------------------------------------------------------------*/

static int
_base_data_4(int        rank_id,
             cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 0;

  _n_g_vtx = 20;

  if (rank_id == 0) {
    _n_vtx = 6;
    n_cells = 2;
    _vtx_range[0] = 0;
    _vtx_range[1] = 5;
  }
  else if (rank_id == 1) {
    _n_vtx = 6;
    n_cells = 2;
    _vtx_range[0] = 5;
    _vtx_range[1] = 10;
  }
  else if (rank_id == 2) {
    _n_vtx = 12;
    n_cells = 5;
    _vtx_range[0] = 10;
    _vtx_range[1] = 14;
  }
  else if (rank_id == 3) {
    _n_vtx = 8;
    n_cells = 3;
    _vtx_range[0] = 14;
    _vtx_range[1] = 20;
  }

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  if (rank_id == 0) {
    _g_vtx_id[0] = 0;
    _g_vtx_id[1] = 1;
    _g_vtx_id[2] = 2;
    _g_vtx_id[3] = 3;
    _g_vtx_id[4] = 4;
    _g_vtx_id[5] = 11;
    cell_vtx[0][0] = 0;
    cell_vtx[0][1] = 1;
    cell_vtx[0][2] = 3;
    cell_vtx[0][3] = 2;
    cell_vtx[1][0] = 2;
    cell_vtx[1][1] = 3;
    cell_vtx[1][2] = 4;
    cell_vtx[1][3] = 5;
  }
  else if (rank_id == 1) {
    _g_vtx_id[0] = 5;
    _g_vtx_id[1] = 6;
    _g_vtx_id[2] = 7;
    _g_vtx_id[3] = 8;
    _g_vtx_id[4] = 9;
    _g_vtx_id[5] = 10;
    cell_vtx[0][0] = 5;
    cell_vtx[0][1] = 0;
    cell_vtx[0][2] = 2;
    cell_vtx[0][3] = 1;
    cell_vtx[1][0] = 1;
    cell_vtx[1][1] = 2;
    cell_vtx[1][2] = 4;
    cell_vtx[1][3] = 3;
  }
  else if (rank_id == 2) {
    _g_vtx_id[0] = 10;
    _g_vtx_id[1] = 11;
    _g_vtx_id[2] = 12;
    _g_vtx_id[3] = 13;
    _g_vtx_id[4] = 1;
    _g_vtx_id[5] = 3;
    _g_vtx_id[6] = 4;
    _g_vtx_id[7] = 6;
    _g_vtx_id[8] = 8;
    _g_vtx_id[9] = 9;
    _g_vtx_id[10] = 14;
    _g_vtx_id[11] = 15;
    cell_vtx[0][0] = 4;
    cell_vtx[0][1] = 0;
    cell_vtx[0][2] = 7;
    cell_vtx[0][3] = 5;
    cell_vtx[1][0] = 5;
    cell_vtx[1][1] = 7;
    cell_vtx[1][2] = 8;
    cell_vtx[1][3] = 6;
    cell_vtx[2][0] = 1;
    cell_vtx[2][1] = 6;
    cell_vtx[2][2] = 2;
    cell_vtx[2][3] = 10;
    cell_vtx[3][0] = 6;
    cell_vtx[3][1] = 8;
    cell_vtx[3][2] = 11;
    cell_vtx[3][3] = 2;
    cell_vtx[4][0] = 8;
    cell_vtx[4][1] = 9;
    cell_vtx[4][2] = 3;
    cell_vtx[4][3] = 11;
  }
  else if (rank_id == 3) {
    _g_vtx_id[0] = 14;
    _g_vtx_id[1] = 15;
    _g_vtx_id[2] = 16;
    _g_vtx_id[3] = 17;
    _g_vtx_id[4] = 18;
    _g_vtx_id[5] = 19;
    _g_vtx_id[6] = 12;
    _g_vtx_id[7] = 13;
    cell_vtx[0][0] = 0;
    cell_vtx[0][1] = 6;
    cell_vtx[0][2] = 3;
    cell_vtx[0][3] = 2;
    cell_vtx[1][0] = 6;
    cell_vtx[1][1] = 1;
    cell_vtx[1][2] = 4;
    cell_vtx[1][3] = 3;
    cell_vtx[2][0] = 1;
    cell_vtx[2][1] = 7;
    cell_vtx[2][2] = 5;
    cell_vtx[2][3] = 4;
  }

  return n_cells;
}

#if defined(REF_4)

static int
_base_data_4_ref(cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 12;

  _n_g_vtx = 20;

  _n_vtx = 20;
  n_cells = 12;
  _vtx_range[0] = 0;
  _vtx_range[1] = 20;

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  for (int i = 0; i < _n_vtx; i++)
    _g_vtx_id[i] = i;

  cell_vtx[0][0] = 0;
  cell_vtx[0][1] = 1;
  cell_vtx[0][2] = 3;
  cell_vtx[0][3] = 2;

  cell_vtx[1][0] = 2;
  cell_vtx[1][1] = 3;
  cell_vtx[1][2] = 4;
  cell_vtx[1][3] = 11;

  cell_vtx[2][0] = 10;
  cell_vtx[2][1] = 5;
  cell_vtx[2][2] = 7;
  cell_vtx[2][3] = 6;

  cell_vtx[3][0] = 6;
  cell_vtx[3][1] = 7;
  cell_vtx[3][2] = 9;
  cell_vtx[3][3] = 8;

  cell_vtx[4][0] = 1;
  cell_vtx[4][1] = 10;
  cell_vtx[4][2] = 6;
  cell_vtx[4][3] = 3;

  cell_vtx[5][0] = 3;
  cell_vtx[5][1] = 6;
  cell_vtx[5][2] = 8;
  cell_vtx[5][3] = 4;

  cell_vtx[6][0] = 11;
  cell_vtx[6][1] = 4;
  cell_vtx[6][2] = 12;
  cell_vtx[6][3] = 14;

  cell_vtx[7][0] = 4;
  cell_vtx[7][1] = 8;
  cell_vtx[7][2] = 15;
  cell_vtx[7][3] = 12;

  cell_vtx[8][0] = 8;
  cell_vtx[8][1] = 9;
  cell_vtx[8][2] = 13;
  cell_vtx[8][3] = 15;

  cell_vtx[9][0] = 14;
  cell_vtx[9][1] = 12;
  cell_vtx[9][2] = 17;
  cell_vtx[9][3] = 16;

  cell_vtx[10][0] = 12;
  cell_vtx[10][1] = 15;
  cell_vtx[10][2] = 18;
  cell_vtx[10][3] = 17;

  cell_vtx[11][0] = 15;
  cell_vtx[11][1] = 13;
  cell_vtx[11][2] = 19;
  cell_vtx[11][3] = 18;

  return n_cells;
}

#endif /* defined(REF_4) */

static int
_base_data_3(int        rank_id,
             cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 0;

  _n_g_vtx = 15;

  if (rank_id == 0) {
    _n_vtx = 6;
    n_cells = 2;
    _vtx_range[0] = 0;
    _vtx_range[1] = 5;
  }
  else if (rank_id == 1) {
    _n_vtx = 6;
    n_cells = 2;
    _vtx_range[0] = 5;
    _vtx_range[1] = 10;
  }
  else if (rank_id == 2) {
    _n_vtx = 12;
    n_cells = 5;
    _vtx_range[0] = 10;
    _vtx_range[1] = 16;
  }

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  if (rank_id == 0) {
    _g_vtx_id[0] = 0;
    _g_vtx_id[1] = 1;
    _g_vtx_id[2] = 2;
    _g_vtx_id[3] = 3;
    _g_vtx_id[4] = 4;
    _g_vtx_id[5] = 11;
    cell_vtx[0][0] = 0;
    cell_vtx[0][1] = 1;
    cell_vtx[0][2] = 3;
    cell_vtx[0][3] = 2;
    cell_vtx[1][0] = 2;
    cell_vtx[1][1] = 3;
    cell_vtx[1][2] = 4;
    cell_vtx[1][3] = 5;
  }
  else if (rank_id == 1) {
    _g_vtx_id[0] = 5;
    _g_vtx_id[1] = 6;
    _g_vtx_id[2] = 7;
    _g_vtx_id[3] = 8;
    _g_vtx_id[4] = 9;
    _g_vtx_id[5] = 10;
    cell_vtx[0][0] = 5;
    cell_vtx[0][1] = 0;
    cell_vtx[0][2] = 2;
    cell_vtx[0][3] = 1;
    cell_vtx[1][0] = 1;
    cell_vtx[1][1] = 2;
    cell_vtx[1][2] = 4;
    cell_vtx[1][3] = 3;
  }
  else if (rank_id == 2) {
    _g_vtx_id[0] = 10;
    _g_vtx_id[1] = 11;
    _g_vtx_id[2] = 12;
    _g_vtx_id[3] = 13;
    _g_vtx_id[4] = 14;
    _g_vtx_id[5] = 15;
    _g_vtx_id[6] = 1;
    _g_vtx_id[7] = 3;
    _g_vtx_id[8] = 4;
    _g_vtx_id[9] = 6;
    _g_vtx_id[10] = 8;
    _g_vtx_id[11] = 9;
    cell_vtx[0][0] = 6;
    cell_vtx[0][1] = 0;
    cell_vtx[0][2] = 9;
    cell_vtx[0][3] = 7;
    cell_vtx[1][0] = 7;
    cell_vtx[1][1] = 9;
    cell_vtx[1][2] = 10;
    cell_vtx[1][3] = 8;
    cell_vtx[2][0] = 1;
    cell_vtx[2][1] = 8;
    cell_vtx[2][2] = 3;
    cell_vtx[2][3] = 2;
    cell_vtx[3][0] = 8;
    cell_vtx[3][1] = 10;
    cell_vtx[3][2] = 4;
    cell_vtx[3][3] = 3;
    cell_vtx[4][0] = 10;
    cell_vtx[4][1] = 11;
    cell_vtx[4][2] = 5;
    cell_vtx[4][3] = 4;
  }

  return n_cells;
}

#if defined(REF_3)

static int
_base_data_3_ref(cs_lnum_t  cell_vtx[][4])
{
  int n_cells = 9;

  _n_g_vtx = 16;
  _n_vtx = _n_g_vtx;
  _vtx_range[0] = 0;
  _vtx_range[1] = 16;

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  for (int i = 0; i < _n_vtx; i++)
    _g_vtx_id[i] = i;

  cell_vtx[0][0] = 0;
  cell_vtx[0][1] = 1;
  cell_vtx[0][2] = 3;
  cell_vtx[0][3] = 2;

  cell_vtx[1][0] = 2;
  cell_vtx[1][1] = 3;
  cell_vtx[1][2] = 4;
  cell_vtx[1][3] = 11;

  cell_vtx[2][0] = 10;
  cell_vtx[2][1] = 5;
  cell_vtx[2][2] = 7;
  cell_vtx[2][3] = 6;

  cell_vtx[3][0] = 6;
  cell_vtx[3][1] = 7;
  cell_vtx[3][2] = 9;
  cell_vtx[3][3] = 8;

  cell_vtx[4][0] = 1;
  cell_vtx[4][1] = 10;
  cell_vtx[4][2] = 6;
  cell_vtx[4][3] = 3;

  cell_vtx[5][0] = 3;
  cell_vtx[5][1] = 6;
  cell_vtx[5][2] = 8;
  cell_vtx[5][3] = 4;

  cell_vtx[6][0] = 11;
  cell_vtx[6][1] = 4;
  cell_vtx[6][2] = 13;
  cell_vtx[6][3] = 12;

  cell_vtx[7][0] = 4;
  cell_vtx[7][1] = 8;
  cell_vtx[7][2] = 14;
  cell_vtx[7][3] = 13;

  cell_vtx[8][0] = 8;
  cell_vtx[8][1] = 9;
  cell_vtx[8][2] = 15;
  cell_vtx[8][3] = 14;

  return n_cells;
}

#endif /* defined(REF_3) */

static int
_base_data_2(int        rank_id,
             cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 0;

  _n_g_vtx = 12;

  if (rank_id == 0) {
    _n_vtx = 8;
    n_cells = 3;
    _vtx_range[0] = 0;
    _vtx_range[1] = 6;
  }
  else if (rank_id == 1) {
    _n_vtx = 8;
    n_cells = 3;
    _vtx_range[0] = 6;
    _vtx_range[1] = 12;
  }

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  if (rank_id == 0) {
    _g_vtx_id[0] = 0;
    _g_vtx_id[1] = 1;
    _g_vtx_id[2] = 2;
    _g_vtx_id[3] = 3;
    _g_vtx_id[4] = 4;
    _g_vtx_id[5] = 5;
    _g_vtx_id[6] = 6;
    _g_vtx_id[7] = 8;
    cell_vtx[0][0] = 0;
    cell_vtx[0][1] = 1;
    cell_vtx[0][2] = 3;
    cell_vtx[0][3] = 2;
    cell_vtx[1][0] = 1;
    cell_vtx[1][1] = 6;
    cell_vtx[1][2] = 7;
    cell_vtx[1][3] = 3;
    cell_vtx[2][0] = 2;
    cell_vtx[2][1] = 3;
    cell_vtx[2][2] = 5;
    cell_vtx[2][3] = 4;
  }
  else if (rank_id == 1) {
    _g_vtx_id[0] = 3;
    _g_vtx_id[1] = 5;
    _g_vtx_id[2] = 6;
    _g_vtx_id[3] = 7;
    _g_vtx_id[4] = 8;
    _g_vtx_id[5] = 9;
    _g_vtx_id[6] = 10;
    _g_vtx_id[7] = 11;
    cell_vtx[0][0] = 0;
    cell_vtx[0][1] = 4;
    cell_vtx[0][2] = 6;
    cell_vtx[0][3] = 1;
    cell_vtx[1][0] = 2;
    cell_vtx[1][1] = 3;
    cell_vtx[1][2] = 5;
    cell_vtx[1][3] = 4;
    cell_vtx[2][0] = 4;
    cell_vtx[2][1] = 5;
    cell_vtx[2][2] = 7;
    cell_vtx[2][3] = 6;
  }

  return n_cells;
}

#if defined(REF_2)

static int
_base_data_2_ref(cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 6;

  _n_g_vtx = 12;
  _n_vtx = _n_g_vtx;
  _vtx_range[0] = 0;
  _vtx_range[1] = 12;

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  for (int i = 0; i < _n_vtx; i++)
    _g_vtx_id[i] = i;

  cell_vtx[0][0] = 0;
  cell_vtx[0][1] = 1;
  cell_vtx[0][2] = 3;
  cell_vtx[0][3] = 2;

  cell_vtx[1][0] = 1;
  cell_vtx[1][1] = 6;
  cell_vtx[1][2] = 8;
  cell_vtx[1][3] = 3;

  cell_vtx[2][0] = 2;
  cell_vtx[2][1] = 3;
  cell_vtx[2][2] = 5;
  cell_vtx[2][3] = 4;

  cell_vtx[3][0] = 3;
  cell_vtx[3][1] = 8;
  cell_vtx[3][2] = 10;
  cell_vtx[3][3] = 5;

  cell_vtx[4][0] = 6;
  cell_vtx[4][1] = 7;
  cell_vtx[4][2] = 9;
  cell_vtx[4][3] = 8;

  cell_vtx[5][0] = 8;
  cell_vtx[5][1] = 9;
  cell_vtx[5][2] = 11;
  cell_vtx[5][3] = 10;

  return n_cells;
}

#endif /* defined(REF_2) */

static int
_base_data_1(cs_lnum_t  cell_vtx[][4])
{
  cs_lnum_t n_cells = 3;

  _n_g_vtx = 8;

  _n_vtx = 8;
  n_cells = 3;

  _vtx_range[0] = 0;
  _vtx_range[1] = 8;

  BFT_MALLOC(_g_vtx_id, _n_vtx, cs_gnum_t);

  _g_vtx_id[0] = 0;
  _g_vtx_id[1] = 1;
  _g_vtx_id[2] = 2;
  _g_vtx_id[3] = 3;
  _g_vtx_id[4] = 4;
  _g_vtx_id[5] = 5;
  _g_vtx_id[6] = 6;
  _g_vtx_id[7] = 7;
  cell_vtx[0][0] = 0;
  cell_vtx[0][1] = 1;
  cell_vtx[0][2] = 5;
  cell_vtx[0][3] = 4;
  cell_vtx[1][0] = 1;
  cell_vtx[1][1] = 2;
  cell_vtx[1][2] = 6;
  cell_vtx[1][3] = 5;
  cell_vtx[2][0] = 2;
  cell_vtx[2][1] = 3;
  cell_vtx[2][2] = 7;
  cell_vtx[2][3] = 6;

  return n_cells;
}

static void
_base_data(int rank_id,
           int n_ranks)
{
  cs_lnum_t n_cells = 0;
  cs_lnum_t cell_vtx[20][4];

  if (n_ranks >= 4)
    n_cells = _base_data_4(rank_id, cell_vtx);
  else if (n_ranks == 3)
    n_cells = _base_data_3(rank_id, cell_vtx);
  else if (n_ranks == 2)
    n_cells = _base_data_2(rank_id, cell_vtx);
  else if (n_ranks == 1) {
    #if defined(REF_4)
    n_cells = _base_data_4_ref(cell_vtx);
    #elif defined(REF_3)
    n_cells = _base_data_3_ref(cell_vtx);
    #elif defined(REF_2)
    n_cells = _base_data_2_ref(cell_vtx);
    #else
    n_cells = _base_data_1(cell_vtx);
    #endif
  }

  /* Compute number of local graph edges */

  _n_edges = n_cells * 4;
  if (_full_connect) _n_edges += n_cells*2;
  if (_symmetric) _n_edges *= 2;

  BFT_MALLOC(_edges, _n_edges, cs_lnum_2_t);

  cs_lnum_t j = 0;
  cs_lnum_t j_step = (_symmetric) ? 2 : 1;
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    _edges[j][0] = cell_vtx[i][0]; _edges[j][1] = cell_vtx[i][1]; j+= j_step;
    _edges[j][0] = cell_vtx[i][1]; _edges[j][1] = cell_vtx[i][2]; j+= j_step;
    _edges[j][0] = cell_vtx[i][2]; _edges[j][1] = cell_vtx[i][3]; j+= j_step;
    _edges[j][0] = cell_vtx[i][3]; _edges[j][1] = cell_vtx[i][0]; j+= j_step;
    if (_full_connect) {
      _edges[j][0] = cell_vtx[i][0]; _edges[j][1] = cell_vtx[i][2]; j+= j_step;
      _edges[j][0] = cell_vtx[i][1]; _edges[j][1] = cell_vtx[i][3]; j+= j_step;
    }
  }
  if (_symmetric) {
    for (cs_lnum_t i = 0; i < _n_edges; i+=2) {
      _edges[i+1][0] = _edges[i][1];
      _edges[i+1][1] = _edges[i][0];
    }
  }
}

/*----------------------------------------------------------------------------
 * build base data based on rank id.
 *----------------------------------------------------------------------------*/

static void
_free_base_data(void)
{
  BFT_FREE(_g_vtx_id);
  BFT_FREE(_edges);
}

/*----------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  /* Internationalization */

#ifdef HAVE_SETLOCALE
  if (!setlocale (LC_ALL,"")) {
#if defined (DEBUG)
     bft_printf("locale not supported by C library"
                " or bad LANG environment variable");
#endif
  }
#endif /* HAVE_SETLOCALE */

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

  _base_data(cs_glob_rank_id, cs_glob_n_ranks);

  cs_gnum_t *g_vtx_num;
  BFT_MALLOC(g_vtx_num, _n_vtx, cs_gnum_t);
  for (cs_lnum_t i = 0; i < _n_vtx; i++)
    g_vtx_num[i] = _g_vtx_id[i] + 1;
  cs_interface_set_t *vtx_ifs
    = cs_interface_set_create(_n_vtx,
                              NULL,
                              g_vtx_num,
                              NULL,
                              0,
                              NULL,
                              NULL,
                              NULL);

  BFT_FREE(g_vtx_num);

  cs_range_set_t *rs
    = cs_range_set_create_from_shared(vtx_ifs,
                                      NULL,
                                      _n_vtx,
                                      _vtx_range,
                                      _g_vtx_id);

  /* Loop on assembler external/internal diagonal */

  for (int id_ie = 0; id_ie < 2; id_ie++) {

    bool sep_diag = (id_ie == 0) ? true : false;

    /* Create a matrix assembler */

    cs_matrix_assembler_t  *ma
      = cs_matrix_assembler_create(_vtx_range, sep_diag);

    cs_matrix_assembler_set_options(ma, 0);

    /* Define connectivities */

    {
      cs_gnum_t g_row_id[20], g_col_id[20];

      /* Diagonal */

      cs_lnum_t j = 0;
      for (cs_lnum_t i = 0; i < _n_vtx; i++) {
        g_row_id[j] = _g_vtx_id[i];
        g_col_id[j] = _g_vtx_id[i];
        j++;
        if (j == 3) {
          cs_matrix_assembler_add_g_ids(ma, j, g_row_id, g_col_id);
          j = 0;
        }
      }
      cs_matrix_assembler_add_g_ids(ma, j, g_row_id, g_col_id);
      j = 0;

      /* Exra-diagonal */
      for (cs_lnum_t i = 0; i < _n_edges; i++) {
        g_row_id[j] = _g_vtx_id[_edges[i][0]];
        g_col_id[j] = _g_vtx_id[_edges[i][1]];
        j++;
        if (j == 3) {
          cs_matrix_assembler_add_g_ids(ma, j, g_row_id, g_col_id);
          j = 0;
        }
      }
      cs_matrix_assembler_add_g_ids(ma, j, g_row_id, g_col_id);
      j = 0;
    }

    /* Now compute structure */

    cs_matrix_assembler_compute(ma);

    /* Assembler is now read for use */

#if 0
    cs_halo_dump(cs_matrix_assembler_get_halo(ma), 1);
#endif

    /* Create associated structures and matrices
       (2 matrices are created simultaneously, to exercice
       the const/shareable aspect of the assembler) */

    cs_matrix_structure_t  *ms_0
      = cs_matrix_structure_create_from_assembler(CS_MATRIX_CSR, ma);
    cs_matrix_structure_t  *ms_1
      = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

    cs_matrix_t  *m_0 = cs_matrix_create(ms_0);
    cs_matrix_t  *m_1 = cs_matrix_create(ms_1);

    /* Now prepare to add values */

    for (int mav_id = 0; mav_id < 2; mav_id++) {

      cs_matrix_assembler_values_t *mav = NULL;

      if (mav_id == 0)
        mav = cs_matrix_assembler_values_init(m_0, 1, 1);
      else
        mav = cs_matrix_assembler_values_init(m_1, 1, 1);

      /* Same ids required as for assembler (at least, no additional ids),
         so loop in a similar manner for safety, but with different
         loop size here (6 instead of 3) */

      cs_gnum_t g_row_id[20], g_col_id[20];
      cs_real_t val[6];
      cs_lnum_t j = 0;

      /* Add terms */

      {
        /* Diagonal */

        for (cs_lnum_t i = 0; i < _n_vtx; i++) {
          g_row_id[j] = _g_vtx_id[i];
          g_col_id[j] = _g_vtx_id[i];
          if (g_row_id[j] < _vtx_range[0] || g_row_id[j] >= _vtx_range[1])
            continue;
          val[j] = cos(g_row_id[j] + 0.1) + sin(g_col_id[j] + 0.1);
          j++;
          if (j == 6) {
            cs_matrix_assembler_values_add_g(mav, j, g_row_id, g_col_id, val);
            j = 0;
          }
        }
        cs_matrix_assembler_values_add_g(mav, j, g_row_id, g_col_id, val);
        j = 0;

        /* Extra-diagonal */

        for (cs_lnum_t i = 0; i < _n_edges; i++) {
          g_row_id[j] = _g_vtx_id[_edges[i][0]];
          g_col_id[j] = _g_vtx_id[_edges[i][1]];
          val[j] = cos(g_row_id[j] + 0.1) + sin(g_col_id[j] + 0.1);
          j++;
          if (j == 6) {
            cs_matrix_assembler_values_add_g(mav, j, g_row_id, g_col_id, val);
            j = 0;
          }
        }
        cs_matrix_assembler_values_add_g(mav, j, g_row_id, g_col_id, val);
        j = 0;
      }

      if (mav_id == 0)
        cs_matrix_assembler_values_done(mav); /* optional */

      cs_matrix_assembler_values_finalize(&mav);

    }

    /* Test SpMV */

    cs_lnum_t n_rows = cs_matrix_get_n_rows(m_0);
    cs_lnum_t n_cols = cs_matrix_get_n_columns(m_0);

    cs_real_t *x, *y_0, *y_1;
    BFT_MALLOC(x, n_cols, cs_real_t);
    BFT_MALLOC(y_0, n_cols, cs_real_t);
    BFT_MALLOC(y_1, n_cols, cs_real_t);
    for (cs_lnum_t i = 0; i < _n_vtx; i++)
      x[i] = (_g_vtx_id[i]+1)*0.5;

    cs_range_set_zero_out_of_range(rs,
                                   CS_REAL_TYPE,
                                   1,
                                   x);

    cs_range_set_gather(rs,
                        CS_REAL_TYPE,
                        1,
                        x,
                        x);

    cs_matrix_vector_multiply(m_0, x, y_0);
    cs_matrix_vector_multiply(m_1, x, y_1);

    bft_printf("\nSpMV pass %d (on range set)\n", id_ie);
    for (cs_lnum_t i = 0; i < n_rows; i++)
      bft_printf("%d: %f %f\n", i, y_0[i], y_1[i]);

    cs_range_set_scatter(rs,
                         CS_REAL_TYPE,
                         1,
                         y_0,
                         y_0);

    cs_range_set_scatter(rs,
                         CS_REAL_TYPE,
                         1,
                         y_1,
                         y_1);

    bft_printf("\nSpMV pass %d (scattered)\n", id_ie);
    for (cs_lnum_t i = 0; i < _n_vtx; i++)
      bft_printf("%d (%d): %f %f\n", i, _g_vtx_id[i], y_0[i], y_1[i]);

    BFT_FREE(x);
    BFT_FREE(y_0);
    BFT_FREE(y_1);

    cs_matrix_release_coefficients(m_0);
    cs_matrix_release_coefficients(m_1);

    cs_matrix_destroy(&m_0);
    cs_matrix_destroy(&m_1);

    cs_matrix_structure_destroy(&ms_0);
    cs_matrix_structure_destroy(&ms_1);

    cs_matrix_assembler_destroy(&ma);
  }

  bft_printf("\n");

  /* Test partition ids on vertices */

  cs_real_t *v0, *v1, *v2;
  BFT_MALLOC(v0, _n_vtx, cs_real_t);
  BFT_MALLOC(v1, _n_vtx, cs_real_t);
  BFT_MALLOC(v2, _n_vtx, cs_real_t);

  for (cs_lnum_t i = 0; i < _n_vtx; i++)
    v0[i] = rs->g_id[i]+1;

  cs_range_set_zero_out_of_range(rs,
                                 CS_REAL_TYPE,
                                 1,
                                 v0);

  for (cs_lnum_t i = 0; i < _n_vtx; i++)
    printf("zero r%d: %d %f\n", cs_glob_rank_id, (int)(rs->g_id[i]), v0[i]);

  bft_printf("\n");

  for (cs_lnum_t i = 0; i < _n_vtx; i++)
    v1[i] = v0[i];

  cs_range_set_gather(rs,
                      CS_REAL_TYPE,
                      1,
                      v0,
                      v2);

  cs_range_set_gather(rs,
                      CS_REAL_TYPE,
                      1,
                      v1,
                      v1);

  for (cs_lnum_t i = 0; i < rs->n_elts[0]; i++)
    printf("gather r%d: %f %f\n", cs_glob_rank_id,
           v1[i], v2[i]);

  bft_printf("\n");

  for (cs_lnum_t i = 0; i < rs->n_elts[1]; i++)
    v0[i] = -2;

  cs_range_set_scatter(rs,
                       CS_REAL_TYPE,
                       1,
                       v2,
                       v0);

  cs_range_set_scatter(rs,
                       CS_REAL_TYPE,
                       1,
                       v1,
                       v1);

  for (cs_lnum_t i = 0; i < rs->n_elts[1]; i++)
    printf("scatter r%d: %d %f %f\n", cs_glob_rank_id,
           (int)(rs->g_id[i]), v0[i], v1[i]);

  BFT_FREE(v2);
  BFT_FREE(v1);
  BFT_FREE(v0);

  cs_range_set_destroy(&rs);

  cs_interface_set_destroy(&vtx_ifs);

  bft_printf("\n");

  /* Finalize */

  _free_base_data();

  bft_mem_end();

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
