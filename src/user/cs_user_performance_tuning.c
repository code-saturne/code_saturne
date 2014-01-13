/*============================================================================
 * Definition of advanced options relative to parallelism.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include "cs_parall.h"
#include "cs_partition.h"
#include "cs_renumber.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define advanced mesh numbering options.
 *----------------------------------------------------------------------------*/

void
cs_user_numbering(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Force the target number of threads for mesh renumbering
     (by default, OMP_NUM_THREADS if OpenMP is enabled, 1 otherwise) */

  if (false)
    cs_renumber_set_n_threads(4);

  /* Set the minimum subset sizes when renumbering for threads. */

  if (false)
    cs_renumber_set_min_subset_size(64,   /* min. interior_subset_size */
                                    64);  /* min. boundary subset_size */

  /* Select renumbering algorithm, among:

     CS_RENUMBER_I_FACES_BLOCK       (no shared cell in block)
     CS_RENUMBER_I_FACES_MULTIPASS   (use multipass face numbering)
     CS_RENUMBER_I_FACES_NONE        (no interior face numbering)
  */

  if (false)
    cs_renumber_set_i_face_algorithm(CS_RENUMBER_I_FACES_MULTIPASS);
}

/*----------------------------------------------------------------------------
 * Define advanced partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_user_partition(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Example:

     Force PT-SCOTCH or SCOTCH for preprocessing partitioning,
     and Hilbert SFC for main partitioning;

     Available algorithms (subject to build with external libraries for
     SCOTCH and METIS) are:

     CS_PARTITION_DEFAULT           Default partitioning, based on stage
     CS_PARTITION_SFC_MORTON_BOX    Morton (Z) curve in bounding box
     CS_PARTITION_SFC_MORTON_CUBE   Morton (Z) curve in bounding cube
     CS_PARTITION_SFC_HILBERT_BOX   Peano-Hilbert curve in bounding box
     CS_PARTITION_SFC_HILBERT_CUBE  Peano-Hilbert curve in bounding cube
     CS_PARTITION_SCOTCH            PT-SCOTCH or SCOTCH
     CS_PARTITION_METIS             ParMETIS or METIS
     CS_PARTITION_BLOCK             Unoptimized (naive) block partitioning */

  if (false) {
    cs_partition_set_algorithm(CS_PARTITION_FOR_PREPROCESS,
                               CS_PARTITION_SCOTCH,
                               1,       /* rank_step */
                               false);  /* ignore periodicity in graph */
    cs_partition_set_algorithm(CS_PARTITION_MAIN,
                               CS_PARTITION_SFC_HILBERT_BOX,
                               1,       /* rank_step */
                               false);  /* ignore periodicity in graph */
  }

  /* Example: set partitioning write to file option.
   *
   * value of write flag:  0: never
   *                       1: for graph-based partitioning only (default)
   *                       2: always */

  if (false)
    cs_partition_set_write_level(0);

  /* Example: force activation/deactivation of initial partitioning
   *          for preprocessing. */

  if (false)
    cs_partition_set_preprocess(false);

  /* Example: define list of extra partitionings to build.
   *
   * Partitionings in this list will be output to file, and may be used for
   * subsequent calculations.
   *
   * When partitioning for both preprocessing and calculation stages, output to
   * file of partioning data or generation of additional partitionings
   * (see \ref cs_partition_add_partitions) will only be done for the
   * second stage. */

  if (false) {
    int  n_extra_partitions = 3;
    int  extra_partitions_list[] = {12, 24, 48};
    cs_partition_add_partitions(n_extra_partitions, extra_partitions_list);
  }
}

/*----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_user_parallel_io(void)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

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

  if (false) {

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

    cs_file_set_mpi_io_positionning(CS_FILE_MPI_INDIVIDUAL_POINTERS);

    MPI_Info_free(&hints);

  }

#endif /* defined(HAVE_MPI_IO) && MPI_VERSION > 1 */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
