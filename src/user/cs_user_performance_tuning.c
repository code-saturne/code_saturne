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

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define advanced mesh numbering options.
 *----------------------------------------------------------------------------*/

void
cs_user_numbering(void)
{
  /* Force the target number of threads for mesh renumbering
     (by default, OMP_NUM_THREADS if OpenMP is enabled, 1 otherwise) */

  if (false)
    cs_renumber_set_n_threads(4);

  /* Set the minimum subset sizes when renumbering for threads. */

  if (false)
    cs_renumber_set_min_subset_size(64,   /* min. interior_subset_size */
                                    64);  /* min. boundary subset_size */

  /* Select renumbering algorithms.

     For cells, available algorithms are:

     CS_RENUMBER_CELLS_SCOTCH_PART     (SCOTCH sub-partitioning, if available)
     CS_RENUMBER_CELLS_SCOTCH_ORDER    (SCOTCH ordering, if available)
     CS_RENUMBER_CELLS_METIS_PART      (METIS sub-partitioning, if available)
     CS_RENUMBER_CELLS_METIS_ORDER     (METIS ordering, if available)
     CS_RENUMBER_CELLS_MORTON          (Morton space filling curve)
     CS_RENUMBER_CELLS_HILBERT         (Hilbert space filling curve)
     CS_RENUMBER_CELLS_NONE            (no renumbering)

     For interior faces, available algorithms are:

     CS_RENUMBER_I_FACES_BLOCK       (no shared cell in block)
     CS_RENUMBER_I_FACES_MULTIPASS   (use multipass face numbering)
     CS_RENUMBER_I_FACES_SIMD        (renumbering for SIMD)
     CS_RENUMBER_I_FACES_NONE        (no interior face numbering)

     Before applying one of those algorithms, interior faces are pre-ordered
     by a lexicographal ordering based on adjacent cells; this ordering
     may be based on the lowest or highest adjacent id first, as defined
     by the CS_RENUMBER_ADJACENT_LOW or CS_RENUMBER_ADJACENT_HIGH value.

     For boundary faces, available algorithms are:

     CS_RENUMBER_B_FACES_THREAD      (renumber for threads)
     CS_RENUMBER_B_FACES_SIMD        (renumbering for SIMD)
     CS_RENUMBER_B_FACES_NONE        (no interior face numbering)
  */

  if (false)
    cs_renumber_set_algorithm
      (false,                           /* halo_adjacent_cells_last */
       false,                           /* halo_adjacent_i_faces_last */
       CS_RENUMBER_ADJACENT_LOW,        /* interior face base ordering  */
       CS_RENUMBER_CELLS_NONE,          /* cells_pre_numbering */
       CS_RENUMBER_CELLS_NONE,          /* cells_numbering */
       CS_RENUMBER_I_FACES_MULTIPASS,   /* interior faces numbering */
       CS_RENUMBER_B_FACES_THREAD);     /* boundary faces numbering */
}

/*----------------------------------------------------------------------------
 * Define advanced partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_user_partition(void)
{
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

/*----------------------------------------------------------------------------
 * Define sparse matrix tuning options.
 *----------------------------------------------------------------------------*/

void
cs_user_matrix_tuning(void)
{
  /* Activate tuning of matrix-vector operations */

  if (false) {

    /* Set tuning runs (defaults) */

    cs_matrix_set_tuning_runs(10,   /* n_min_products */
                              0.5); /* t_measure */

    /* Activate tuning for selected matrix fill types. */

    cs_matrix_set_tuning(CS_MATRIX_SCALAR, 1);

    cs_matrix_set_tuning(CS_MATRIX_SCALAR_SYM, 1);

    /* Force variant for selected types */

    cs_matrix_variant_t *mv
      = cs_matrix_variant_create(CS_MATRIX_MSR,
                                 cs_glob_mesh->i_face_numbering);
    cs_matrix_variant_set_func(mv,
                               cs_glob_mesh->i_face_numbering,
                               CS_MATRIX_33_BLOCK_D,
                               2,
                               "default");

    cs_matrix_set_variant(CS_MATRIX_33_BLOCK_D, mv);

    cs_matrix_variant_destroy(&mv);

    /* Also allow tuning for multigrid for all expected levels
     * (we rarely have more than 10 or 11 levels except for huge meshes). */

    cs_grid_set_matrix_tuning(CS_MATRIX_SCALAR_SYM, 12);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
