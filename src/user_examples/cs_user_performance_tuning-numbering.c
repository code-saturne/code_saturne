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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_performance_tuning-numbering.c
 *
 * \brief Mesh numbering example.
 *
 * See \subpage cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced mesh numbering options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_numbering(void)
{
  /*! [performance_tuning_numbering] */

  /* Force the target number of threads for mesh renumbering
     (by default, OMP_NUM_THREADS if OpenMP is enabled, 1 otherwise) */

  cs_renumber_set_n_threads(4);

  /* Set the minimum subset sizes when renumbering for threads. */

  cs_renumber_set_min_subset_size(64,   /* min. interior_subset_size */
                                  64);  /* min. boundary subset_size */

  /* Select renumbering algorithms */

  cs_renumber_set_algorithm
    (false,                           /* halo_adjacent_cells_last */
     false,                           /* halo_adjacent_i_faces_last */
     CS_RENUMBER_ADJACENT_LOW,        /* interior face base ordering  */
     CS_RENUMBER_CELLS_NONE,          /* cells_pre_numbering */
     CS_RENUMBER_CELLS_NONE,          /* cells_numbering */
     CS_RENUMBER_I_FACES_MULTIPASS,   /* interior faces numbering */
     CS_RENUMBER_B_FACES_THREAD,      /* boundary faces numbering */
     CS_RENUMBER_VERTICES_NONE);      /* vertices numbering */

  /*! [performance_tuning_numbering] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
