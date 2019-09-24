/*============================================================================
 * Definition of advanced options relative to parallelism.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_file.h"
#include "cs_grid.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_parall.h"
#include "cs_partition.h"
#include "cs_renumber.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_performance_tuning-parallel-io.c
 *
 * \brief Parallel Input/Output example.
 *
 * See \subpage cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define parallel IO settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parallel_io(void)
{
  /*! [perfomance_tuning_parallel_io] */

#if defined(HAVE_MPI_IO) && MPI_VERSION > 1

  /* Example fine-tune parallel IO settings.

     Available distributed block access methods
     (subject to build with MPI IO) are:

     CS_FILE_STDIO_SERIAL        Serial standard C IO
                                 (funnelled through rank 0 in parallel)
     CS_FILE_STDIO_PARALLEL      Per-process standard C IO
     CS_FILE_MPI_INDEPENDENT     Non-collective MPI-IO
                                 with independent file open and close
     CS_FILE_MPI_NON_COLLECTIVE  Non-collective MPI-IO
                                 with collective file open and close
     CS_FILE_MPI_COLLECTIVE      Collective MPI-IO
  */

  int block_rank_step = 8;
  int block_min_size = 1024*1024*8;

  MPI_Info hints = MPI_INFO_NULL;
  cs_file_access_t  method = CS_FILE_MPI_COLLECTIVE;

  /* Set MPI IO hints

     (see MPI-IO or your filesystem documentation;
     examples here may have no effect, improve, or degrade performance)

     For LUSTRE filesystems, many articles in the literature seem
     to recommend adjusting striping to improve performance.

     If using ROMIO, useful hints for collective buffering and data-sieving
     may take values: "enable", "disable", "automatic".
  */

  MPI_Info_create(&hints);

  MPI_Info_set(hints, "striping_factor", "8");
  MPI_Info_set(hints, "striping_unit", "8388608");

  MPI_Info_set(hints, "romio_cb_read", "automatic");
  MPI_Info_set(hints, "romio_cb_write", "automatic");
  MPI_Info_set(hints, "romio_ds_read", "automatic");
  MPI_Info_set(hints, "romio_ds_write", "automatic");

  /* Set default file acces methods and communicator stride */

  cs_file_set_default_access(CS_FILE_MODE_WRITE, method, hints);

  MPI_Info_set(hints, "collective_buffering", "true");
  MPI_Info_set(hints, "access_style", "read_once");

  cs_file_set_default_access(CS_FILE_MODE_READ, method, hints);

  cs_file_set_default_comm(block_rank_step, block_min_size, cs_glob_mpi_comm);

  cs_file_set_mpi_io_positioning(CS_FILE_MPI_INDIVIDUAL_POINTERS);

  MPI_Info_free(&hints);

#endif /* defined(HAVE_MPI_IO) && MPI_VERSION > 1 */

  /*! [perfomance_tuning_parallel_io] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
