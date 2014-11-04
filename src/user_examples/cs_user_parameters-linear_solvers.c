/*============================================================================
 * User subroutines for input of calculation parameters.
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
#include <string.h>
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism alowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{
  /* Available native iterative linear solvers are:
   *
   *  CS_SLES_PCG        (preconditionned conjugate gradient)
   *  CS_SLES_JACOBI     (Jacobi)
   *  CS_SLES_BICGSTAB   (Bi-conjugate gradient stabilized)
   *  CS_SLES_BICGSTAB2  (BiCGStab2)
   *  CS_SLES_GMRES      (generalized minimal residual)
   *
   *  The multigrid solver uses the conjugate gradient as a smoother
   *  and coarse solver by default, but this behavior may be modified. */

  /* Example: use multigrid for wall distance computation */
  /*------------------------------------------------------*/

  /*! [sles_wall_dist] */
  cs_multigrid_define(-1, "wall_distance");
  /*! [sles_wall_dist] */

  /* Example: use BiCGStab2 for user variable (named user_1) */
  /*---------------------------------------------------------*/

  /*! [sles_user_1] */
  cs_field_t *cvar_user_1 = cs_field_by_name_try("user_1");
  if (cvar_user_1 != NULL) {
    cs_sles_it_define(cvar_user_1->id,
                      NULL,   /* name passed is NULL if field_id > -1 */
                      CS_SLES_BICGSTAB2,
                      1,      /*  polynomial precond. degree (default 0) */
                      10000); /* n_max_iter */
  }
  /*! [sles_user_1] */

  /* Example: increase verbosity parameters for pressure */

  /*! [sles_verbosity_1] */

  cs_sles_t *sles_p = cs_sles_find_or_add(CS_F_(p)->id, NULL);
  cs_sles_set_verbosity(sles_p, 4);

  /*! [sles_verbosity_1] */

  /* Example: change multigrid parameters for pressure */

  /*! [sles_mgp_1] */
  cs_multigrid_t *mg = cs_multigrid_define(CS_F_(p)->id, NULL);

  cs_multigrid_set_coarsening_options(mg,
                                      3,    /* aggregation_limit (default 3) */
                                      0,    /* coarsening_type (default 0) */
                                      10,   /* n_max_levels (default 25) */
                                      30,   /* min_g_cells (default 30) */
                                      0.95, /* P0P1 relaxation (default 0.95) */
                                      20);  /* postprocessing (default 0) */

  cs_multigrid_set_solver_options
    (mg,
     CS_SLES_JACOBI, /* descent smoother type (default: CS_SLES_PGC) */
     CS_SLES_JACOBI, /* ascent smoother type (default: CS_SLES_PGC) */
     CS_SLES_PCG,    /* coarse solver type (default: CS_SLES_PGC) */
     50,             /* n max cycles (default 100) */
     5,              /* n max iter for descent (default 10) */
     5,              /* n max iter for asscent (default 10) */
     1000,           /* n max iter coarse solver (default 10) */
     0,              /* polynomial precond. degree descent (default 0) */
     0,              /* polynomial precond. degree ascent (default 0) */
     1,              /* polynomial precond. degree coarse (default 0) */
     -1.0,           /* precision multiplier descent (< 0 forces max iters) */
     -1.0,           /* precision multiplier ascent (< 0 forces max iters) */
     0.1);           /* requested precision multiplier coarse (default 1) */
  /*! [sles_mgp_1] */

  /* Set parallel grid merging options for all multigrid solvers */

  /*! [sles_mg_parall] */
  cs_grid_set_merge_options(4,   /* # of ranks merged at a time */
                            300, /* mean # of cells under which we merge */
                            500, /* global # of cells under which we merge */
                            1);  /* # of ranks under which we do not merge */
  /*! [sles_mg_parall] */

  /* Set a non-default linear solver for DOM radiation.
     The solver must be set for each direction; here, we assume
     a quadrature with 32 directions is used */

  /*! [sles_rad_dom_1] */
  for (int i = 0; i < 32; i++) {
    char name[16];
    sprintf(name, "radiation_%03d", i+1);
    cs_sles_it_define(-1,
                      name,
                      CS_SLES_JACOBI,
                      0,      /* poly_degree */
                      1000);  /* n_max_iter */

  }
  /*! [sles_rad_dom_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
