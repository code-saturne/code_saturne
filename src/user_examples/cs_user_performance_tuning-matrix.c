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
 * \file cs_user_performance_tuning-matrix.c
 *
 * \brief Matrix tuning example.
 *
 * See \subpage cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse matrix tuning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_matrix_tuning(void)
{
  /*! [performance_tuning_matrix] */

  /* Activate tuning of matrix-vector operations */

  /* Set tuning runs (defaults) */

  cs_matrix_set_tuning_runs(20,   /* n_min_products */
                            0.5); /* t_measure */

  /* Force default for selected types */

  cs_matrix_default_set_type(CS_MATRIX_BLOCK_D, CS_MATRIX_MSR);

  /* Also allow tuning for multigrid for all expected levels
   * (we rarely have more than 10 or 11 levels except for huge meshes). */

  cs_grid_set_matrix_tuning(CS_MATRIX_SCALAR_SYM, 12);

  /*! [performance_tuning_matrix] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
