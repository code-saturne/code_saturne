/*============================================================================
 * Data checking for the 1D thermal wall module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_restart.h"
#include "cs_1d_wall_thermal.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_1d_wall_thermal_check.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Data checking for the 1D thermal wall module.
 *
 * \param[in]   iappel   Call number:
 *                       - 1: first call during the initialization (called once)
 *                       Checking the number of cells nfpt1d and isuit1.
 *                       - 2: second call during the initialization (called once)
 *                       Checking ifpt1d, nppt1d, eppt1d and rgpt1d arrays.
 *                       - 3: called at each time step
 *                       Checking iclt1d, xlmbt1, rcpt1d and dtpt1d arrays.
 * \param[in]   isuit1   Rereading of the restart file:
 *                       - 0: No rereading
 *                            (meshing and wall temperature reinitialization)
 *                       - 1: Rereading of the restart file for the 1-Dimension
 *                            thermal module
 *                       - isuite: Rereading only if the computational fluid
 *                            dynamic is a continuation of the computation.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_check(int iappel,
                         int isuit1)
{
  cs_lnum_t ii, ifac;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_lnum_t nfpt1d = cs_glob_1d_wall_thermal->nfpt1d;

  if (iappel == 1) {
    if (nfpt1d < 0 || nfpt1d > n_b_faces)  {
      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                 "@    ========\n"
                 "@    1D-WALL THERMAL MODULE\n"
                 "@\n"
                 "@    NFPT1D MUST BE POSITIVE AND LOWER THAN NFABOR\n"
                 "@    ONE HAS HERE\n"
                 "@       NFABOR = %d\n"
                 "@       NFPT1D = %d\n"
                 "@\n"
                 "@  The calculation will not run.\n"
                 "@\n"
                 "@  Verify uspt1d.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n", n_b_faces, nfpt1d);
      cs_exit(EXIT_FAILURE);
    }
    if (isuit1 != 0 && isuit1 != 1) {
      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                 "@    ========\n"
                 "@    1D-WALL THERMAL MODULE\n"
                 "@\n"
                 "@    ISUIT1 MUST BE AN INTEGER EQUAL TO 0 OR 1\n"
                 "@    ITS VALUE IS %d\n"
                 "@\n"
                 "@  The calculation will not run.\n"
                 "@\n"
                 "@  Verify uspt1d.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n", isuit1);
      cs_exit(EXIT_FAILURE);
    }

  } else if (iappel == 2) {
      for (ii = 0 ; ii < nfpt1d ; ii++) {
        cs_lnum_t ifpt1d = cs_glob_1d_wall_thermal->ifpt1d[ii] - 1;
        if (ifpt1d < 0 || ifpt1d > n_b_faces) {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY IFPT1D MUST GIVE A BOUNDARY FACE NUMBER\n"
                     "@    ONE HAS HERE\n"
                     "@       NFABOR = %d\n"
                     "@       IFPT1D(%d) = %d\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", n_b_faces, ii, ifpt1d);
          cs_exit(EXIT_FAILURE);
        }
      }

      for (ii = 0 ; ii < nfpt1d ; ii++) {
        ifac = cs_glob_1d_wall_thermal->ifpt1d[ii] - 1;
        if (cs_glob_1d_wall_thermal->local_models[ii].nppt1d <= 0)  {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY NPPT1D MUST GIVE A POSITIVE INTEGER\n"
                     "@    ONE HAS HERE\n"
                     "@       NPPT1D(%d) = %d\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@", ii, cs_glob_1d_wall_thermal->local_models[ii].nppt1d);
          cs_exit(EXIT_FAILURE);
        }
        if (cs_glob_1d_wall_thermal->local_models[ii].eppt1d <= 0.) {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY EPPT1D MUST GIVE A POSITIVE REAL\n"
                     "@    ONE HAS HERE\n"
                     "@       EPPT1D(%d) = %14.5e\n"
                     "@       (BOUNDARY FACE NUMBER %d)\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", ii, cs_glob_1d_wall_thermal->local_models[ii].eppt1d, ifac);
          cs_exit(EXIT_FAILURE);
        }
        if (cs_glob_1d_wall_thermal->local_models[ii].rgpt1d <= 0.) {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY RGPT1D MUST GIVE A POSITIVE REAL\n"
                     "@    ONE HAS HERE\n"
                     "@       RGPT1D(%d) = %14.5e\n"
                     "@       (BOUNDARY FACE NUMBER %d)\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", ii, cs_glob_1d_wall_thermal->local_models[ii].rgpt1d, ifac);
          cs_exit(EXIT_FAILURE);
        }
      }

  } else if (iappel == 3) {
    for (ii = 0 ; ii < nfpt1d ; ii++) {
      ifac = cs_glob_1d_wall_thermal->ifpt1d[ii] - 1;
      int iclt1d = cs_glob_1d_wall_thermal->local_models[ii].iclt1d;
      if (iclt1d != 1 && iclt1d != 3)  {
        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                   "@    ========\n"
                   "@    1D-WALL THERMAL MODULE\n"
                   "@\n"
                   "@    THE ARRAY ICLT1D CAN ONLY TAKE THE VALUES 1 OR 3\n"
                   "@    ONE HAS HERE\n"
                   "@       ICLT1D(%d) = %d\n"
                   "@       (BOUNDARY FACE NUMBER %d)\n"
                   "@\n"
                   "@  The calculation will not run.\n"
                   "@\n"
                   "@  Verify uspt1d.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n", ii, iclt1d, ifac);
        cs_exit(EXIT_FAILURE);
      }
      if (cs_glob_1d_wall_thermal->local_models[ii].xlmbt1 <= 0.)  {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY XLMBT1 MUST GIVE A POSITIVE REAL\n"
                     "@    ONE HAS HERE\n"
                     "@       XLMBT1(%d) = %14.5e\n"
                     "@       (BOUNDARY FACE NUMBER %d)\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", ii, cs_glob_1d_wall_thermal->local_models[ii].xlmbt1, ifac);
        cs_exit(EXIT_FAILURE);
      }
      if (cs_glob_1d_wall_thermal->local_models[ii].rcpt1d <= 0.)  {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY RCPT1D MUST GIVE A POSITIVE REAL\n"
                     "@    ONE HAS HERE\n"
                     "@       RCPT1D(%d) = %14.5e\n"
                     "@       (BOUNDARY FACE NUMBER %d)\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", ii, cs_glob_1d_wall_thermal->local_models[ii].rcpt1d, ifac);
        cs_exit(EXIT_FAILURE);
      }
      if (cs_glob_1d_wall_thermal->local_models[ii].dtpt1d <= 0.)  {
          bft_printf("@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n"
                     "@ @@ WARNING: ABORT DURING THE DATA SPECIFICATION\n"
                     "@    ========\n"
                     "@    1D-WALL THERMAL MODULE\n"
                     "@\n"
                     "@    THE ARRAY DTPT1D MUST GIVE A POSITIVE REAL\n"
                     "@    ONE HAS HERE\n"
                     "@       DTPT1D(%d) = %14.5e\n"
                     "@       (BOUNDARY FACE NUMBER %d)\n"
                     "@\n"
                     "@  The calculation will not run.\n"
                     "@\n"
                     "@  Verify uspt1d.\n"
                     "@\n"
                     "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                     "@\n", ii, cs_glob_1d_wall_thermal->local_models[ii].dtpt1d, ifac);
        cs_exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
