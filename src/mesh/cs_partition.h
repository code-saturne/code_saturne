#ifndef __CS_PARTITION_H__
#define __CS_PARTITION_H__

/*============================================================================
 * Define cs_mesh_t fields from cs_mesh_builder_t fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_base.h"
#include "cs_log.h"

#include "cs_mesh.h"
#include "cs_mesh_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Partitioning stage
 *
 * Partitioning is always done just after reading the mesh, unless a
 * partitioning input file is available, in which case the partitioning
 * read replaces this stage.
 *
 * When a mesh modification implying a change of cell connectivity graph
 * is expected, the mesh may be re-partitioned after the pre-processing
 * stage, prior to calculation. By default, re-partitioning is only done
 * if the partitioning algorithm chosen for that stage is expected to
 * produce different results due to the connectivity change. This is
 * the case for graph-based algorithms such as those of METIS or SCOTCH,
 * when mesh joining is defined, or additional periodic matching is defined
 * (and the algorithm is not configured to ignore periodicity information).
 *
 * There are thus two possible partitioning stages:
 *
 * - CS_PARTITION_FOR_PREPROCESS, which is optional, and occurs
 *   just  after reading the mesh.
 * - CS_PARTITION_MAIN, which occurs just after reading the mesh if
 *   it is the only stage,, or after mesh preprocessing (and before
 *   computation), if the partitioning for preprocessing stage is
 *   activated.
 *
 * The number of partitioning stages is determined automatically based on
 * information provided through cs_partition_set_preprocess_hints(),
 * but re-partitioning may also be forced or inhibited using the
 * cs_partition_set_preprocess() function.
 */

typedef enum {

  CS_PARTITION_FOR_PREPROCESS,  /* Partitioning for preprocessing stage */
  CS_PARTITION_MAIN             /* Partitioning for computation stage */

} cs_partition_stage_t;


/* Partitioning algorithm type
 *
 * If the default algorithm is selected, the choice will be based on the
 * following priority, depending on available libraries:
 * -  Pt-Scotch (or Scotch if partitioning on one rank);
 * -  ParMETIS (or METIS if partitioning on one rank);
 * -  Morton space-filling curve (in bounding box)
 *
 * If both partitioning stages are active, the default for the preprocessing
 * stage will be based on the Morton space-filling curve (in bounding box),
 * as this should be cheaper, and the initial cell connectivity graph
 * is usually expected to be modified during preprocessing.
 */

typedef enum {

  CS_PARTITION_DEFAULT,           /* Default partitioning (based on stage) */
  CS_PARTITION_SFC_MORTON_BOX,    /* Morton (Z) curve in bounding box */
  CS_PARTITION_SFC_MORTON_CUBE,   /* Morton (Z) curve in bounding cube */
  CS_PARTITION_SFC_HILBERT_BOX,   /* Peano-Hilbert curve in bounding box */
  CS_PARTITION_SFC_HILBERT_CUBE,  /* Peano-Hilbert curve in bounding cube */
  CS_PARTITION_SCOTCH,            /* PT-SCOTCH or SCOTCH */
  CS_PARTITION_METIS,             /* ParMETIS or METIS */
  CS_PARTITION_BLOCK              /* Unoptimized (naive) block partitioning */

} cs_partition_algorithm_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print information on external libraries
 *
 * parameters:
 *   log_type  <--  log type
 *----------------------------------------------------------------------------*/

void
cs_partition_external_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------
 * Set algorithm for domain partitioning for a given partitioning stage.
 *
 * parameters:
 *   stage        <-- associated partitioning stage
 *   algorithm    <-- partitioning algorithm choice
 *   rank_step    <-- if > 1, partitioning done on at most
 *                    n_ranks / rank_step processes
 *                    (for graph-based partitioning only)
 *   ignore_perio <-- if true, ignore periodicity information when present
 *                    when present (for graph-based
 *                    (for graph-based partitioning only)
 *----------------------------------------------------------------------------*/

void
cs_partition_set_algorithm(cs_partition_stage_t      stage,
                           cs_partition_algorithm_t  algorithm,
                           int                       rank_step,
                           bool                      ignore_perio);

/*----------------------------------------------------------------------------
 * Set partitioning write to file option.
 *
 * Partitioning information for subsequent calculations is written to file
 * after the last partitioning stage depending on the output level.
 *
 * Note that partitioning information for additional partitionings is
 * always written to file, regardless of this option.
 *
 * parameters:
 *   write_flag <-- option to save partitioning information:
 *                  0: never
 *                  1: for graph-based partitioning only (default)
 *                  2: always
 *----------------------------------------------------------------------------*/

void
cs_partition_set_write_level(int  write_flag);

/*----------------------------------------------------------------------------
 * Define hints indicating if initial partitioning fo a preprocessing
 * stage is required.
 *
 * parameters:
 *   join          <-- true if a mesh joining operation is planned
 *   join_periodic <-- true if a mesh periodic matching operation is planned
 *----------------------------------------------------------------------------*/

void
cs_partition_set_preprocess_hints(bool  join,
                                  bool  join_periodic);

/*----------------------------------------------------------------------------
 * Activate or deactivate initial partitioning for preprocessing.
 *
 * parameters:
 *   active <-- true to activate pre-partitiong for the preprocessing
 *              stage, false to de-activate it
 *----------------------------------------------------------------------------*/

void
cs_partition_set_preprocess(bool  active);

/*----------------------------------------------------------------------------
 * Indicate if initial partitioning for preprocessing is required.
 *
 * returns:
 *   true if initial partitioning for preprocessing is active,
 *   false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_partition_get_preprocess(void);

/*----------------------------------------------------------------------------
 * Define list of extra partitionings to build.
 *
 * Partitionings in this list will be output to file, and may be used for
 * subsequent calculations.
 *
 * When partitioning for both preprocessing and calculation stages, output to
 * file of partioning data or generation of additional partitionings
 * (see \ref cs_partition_add_partitions) will only be done for the
 * second stage.
 *
 * parameters:
 *   n_extra_partitions    <-- number of extra partitionings to compute
 *   extra_partitions_list <-- list of extra partitions to compute
 *----------------------------------------------------------------------------*/

void
cs_partition_add_partitions(int  n_extra_partitions,
                            int  extra_partitions_list[]);

/*----------------------------------------------------------------------------
 * Compute partitioning for a given mesh.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-> pointer to mesh builder structure
 *   stage        <-- associated partitioning stage
 *----------------------------------------------------------------------------*/

void
cs_partition(cs_mesh_t             *mesh,
             cs_mesh_builder_t     *mesh_builder,
             cs_partition_stage_t   stage);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARTITION_H__ */
