/*============================================================================
 * Data Entry of the 1D wall thermal module.
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_1d_wall_thermal.c
 *
 * \brief Data Entry of the 1D wall thermal module.
 *
 * \param[in]   iappel   Call number:
 *                       - 1: first call during the initialization (called once)
 *                       Setting the number of cells nfpt1d.
 *                       - 2: second call during the initialization (called once)
 *                       Filling ifpt1d, nppt1d, eppt1d and rgpt1d arrays.
 *                       - 3: called at each time step
 *                       Filling iclt1d, xlmbt1, rcpt1d and dtpt1d arrays.

 * \param[in]   isuit1   Rereading of the restart file:
 *                       - 0: No rereading
 *                            (meshing and wall temperature reinitialization)
 *                       - 1: Rereading of the restart file for the 1-Dimension
 *                            thermal module
 *                       - isuite: Rereading only if the computational fluid
 *                            dynamic is a continuation of the computation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

void
cs_user_1d_wall_thermal(int iappel,
                        int isuit1)
{
  /*! [loc_var_dec] */
  int izone, ifbt1d;
  cs_lnum_t nlelt;
  cs_lnum_t *lstelt;
  cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_real_3_t *cdgfbo = (cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  cs_1d_wall_thermal_t *wall_thermal = cs_get_glob_1d_wall_thermal();

/*! [loc_var_dec] */

/*! [allocate] */

  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

/*! [allocate] */

/*! [restart] */

  /*----------------------------------------------------------------------------*
   * Rereading of the restart file:
   *----------------------------------

   *    isuit1 = 0        --> No rereading
   *                          (meshing and wall temperature reinitialization)
   *    isuit1 = 1        --> Rereading of the restart file for the 1-Dimension
   *                          thermal module
   *    isuit1 = isuite   --> Rereading only if the computational fluid dynamic is
   *                           a continuation of the computation.
   *    The initialization of isuit1 is mandatory.
   *----------------------------------------------------------------------------*/

  isuit1 = cs_restart_present();
  izone = 0;
  ifbt1d = 0;

/*! [restart] */

/*! [iappel_12] */

  if (iappel == 1 || iappel == 2) {
    /*----------------------------------------------------------------------------*
     * Faces determining with the 1-D thermal module:
     *----------------------------------------------
     *
     *     nfpt1d    : Total number of faces with the 1D thermal module
     *     ifpt1d[ii]: Number of the [ii]th face with the 1D thermal module

     * Remarks:
     *--------
     *     During the rereading of the restart file, nfpt1d and ifpt1d are
     *     compared with the other values from the restart file being the result of
     *     the start or restarting computation.
     *
     *     A total similarity is required to continue with the previous computation.
     *     Regarding the test case on ifpt1d, it is necessary that the array be
     *     arranged in increasing order (ifpt1[jj] > ifpt1d[ii] if jj > ii).
     *
     *     If it is impossible, contact the developer team to deactivate this test.
     *----------------------------------------------------------------------------*/

    /* Get the list of boundary faces that will be coupled */

    cs_selector_get_b_face_num_list("2 or 3 or 5 or 6 or 7 or 8 or 9 or 10",
                                    &nlelt, lstelt);

    izone++;

    /* Fill the ifpt1d array, and compute nfpt1d */

    for (cs_lnum_t ilelt = 0 ; ilelt < nlelt ; ilelt++) {
      cs_lnum_t ifac = lstelt[ilelt];
      wall_thermal->izft1d[ifac-1] = izone;
      if (iappel == 2) wall_thermal->ifpt1d[ifbt1d] = ifac;
      ifbt1d++;
    }
  }

  if (iappel == 1) {
    wall_thermal->nfpt1d = ifbt1d;
  }

/*! [iappel_12] */

/*! [iappel_2] */

  /*----------------------------------------------------------------------------*
   * Parameters padding of the mesh and initialization:
   *--------------------------------------------------
   *
   *     (Only one pass during the beginning of the computation)

   *     locals_models[ii].nppt1d: number of discretized points associated
   *                 to the (ii)th face with the 1-D thermal module.
   *     locals_models[ii].eppt1d: wall thickness associated to the (ii)th face
   *                 with the 1-D thermal module.
   *     locals_models[ii].rgpt1d: geometric progression ratio of the meshing refinement
   *                 associated to the (ii)th face with the 1-D thermal module.
   *                 (with : rgpt1d > 1 => small meshes on the fluid side)
   *     locals_models[ii].tppt1d: wall temperature initialization associated to the
   *                 (ii)th face with the 1-D thermal module.

   * Remarks:
   *--------
   *     During the rereading of the restart file for the 1-D thermal module,
   *     the tppt1d variable is not used.
   *
   *     The nfpt1d, eppt1d and rgpt1d variables are compared to the previous
   *     values being the result of the restart file.
   *
   *     An exact similarity is necessary to continue with the previous computation.
   *----------------------------------------------------------------------------*/

  if (iappel == 2) {
    for (cs_lnum_t ii = 0 ; ii < wall_thermal->nfpt1d ; ii++) {
      wall_thermal->local_models[ii].nppt1d = 8;
      wall_thermal->local_models[ii].eppt1d = 0.01144;
      wall_thermal->local_models[ii].rgpt1d = 1.;
      wall_thermal->tppt1d[ii] = cs_glob_fluid_properties->t0;
    }
  }

/*! [iappel_2] */

/*! [iappel_3] */

  /*----------------------------------------------------------------------------*
   * Padding of the wall exterior boundary conditions:
   *-------------------------------------------------
   *
   *    local_models[ii].iclt1d: boundary condition type
   *     ----------
   *    local_models[ii].iclt1d = 1: dirichlet condition,
   *                                  with exchange coefficient
   *    local_models[ii].iclt1d = 3: flux condition
   *
   *    local_models[ii].tept1d: exterior temperature
   *    local_models[ii].hept1d: exterior exchange coefficient
   *    local_models[ii].fept1d: flux applied to the exterior (flux<0 = coming flux)
   *    local_models[ii].xlmt1d: lambda wall conductivity coefficient (W/m/C)
   *    local_models[ii].rcpt1d: wall coefficient rho*Cp (J/m3/C)
   *    local_models[ii].dtpt1d: time step resolution of the thermal equation to the
   *                             (ii)th boundary face with the 1-D thermal module (s)
   *----------------------------------------------------------------------------*/

  if (iappel == 3) {
    for (cs_lnum_t ii = 0 ; ii < wall_thermal->nfpt1d ; ii++) {
      wall_thermal->local_models[ii].iclt1d = 1;
      wall_thermal->local_models[ii].tept1d = cs_glob_fluid_properties->t0;
      /* Physical parameters */
      cs_lnum_t face_id = wall_thermal->ifpt1d[ii] - 1;
      /* Floor: plate */
      if (cdgfbo[face_id][2] <= 1.e-3) {
        wall_thermal->local_models[ii].xlmbt1 = 0.16;
        wall_thermal->local_models[ii].rcpt1d = 790.*900.;
      /* Wall and ceiling: marinite */
      } else {
        wall_thermal->local_models[ii].xlmbt1 = 0.11;
        wall_thermal->local_models[ii].rcpt1d = 670.*778.;
      }
      cs_lnum_t c_id = b_face_cells[face_id];
      wall_thermal->local_models[ii].dtpt1d = CS_F_(dt)->val[c_id];
    }
  }

/*! [iappel_3] */

  /*----------------------------------------------------------------------------*
   * End of the cs_user_1d_wall_thermal function
   *----------------------------------------------------------------------------*/

/*! [deallocate] */

  BFT_FREE(lstelt);

/*! [deallocate] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
