/*============================================================================
 * Definition of advanced options relative to parallelism.
 *============================================================================*/

/* VERS */

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_performance_tuning-partition.c
 *
 * \brief Partitioning option setting example.
 *
 * See \ref cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced partitioning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_partition(void)
{
  /*! [performance_tuning_partition_1] */
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

    cs_partition_set_algorithm(CS_PARTITION_FOR_PREPROCESS,
                               CS_PARTITION_SCOTCH,
                               1,       /* rank_step */
                               false);  /* ignore periodicity in graph */
    cs_partition_set_algorithm(CS_PARTITION_MAIN,
                               CS_PARTITION_SFC_HILBERT_BOX,
                               1,       /* rank_step */
                               false);  /* ignore periodicity in graph */

  }
  /*! [performance_tuning_partition_1] */

  /*! [performance_tuning_partition_2] */
  {
    /* Example: set partitioning write to file option.
     *
     * value of write flag:  0: never
     *                       1: for graph-based partitioning only (default)
     *                       2: always */

    cs_partition_set_write_level(0);

  }
  /*! [performance_tuning_partition_2] */


  /*! [performance_tuning_partition_3] */
  {
    /* Example: force activation/deactivation of initial partitioning
     *          for preprocessing. */

    cs_partition_set_preprocess(false);
  }
  /*! [performance_tuning_partition_3] */


  /*! [performance_tuning_partition_4] */
  {
    /* Example: define list of extra partitionings to build.
     *
     * Partitionings in this list will be output to file, and may be used for
     * subsequent calculations.
     *
     * When partitioning for both preprocessing and calculation stages, output to
     * file of partioning data or generation of additional partitionings
     * (see \ref cs_partition_add_partitions) will only be done for the
     * second stage. */

    int  n_extra_partitions = 3;
    int  extra_partitions_list[] = {12, 24, 48};
    cs_partition_add_partitions(n_extra_partitions, extra_partitions_list);
  }
  /*! [performance_tuning_partition_4] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
