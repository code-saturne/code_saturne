/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between Code_Saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * \file cs_user_mesh-quality.c
 *
 * \brief Mesh quality example
 *
 * See \subpage cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set options for cutting of warped faces
 *
 * \param[in,out] mesh pointer to cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_warping(void)
{
  /*! [mesh_warping] */

  double max_warp_angle = 3; /* bounded between 0 and 90 degrees */
  int postprocess = 0;

  cs_mesh_warping_set_defaults(max_warp_angle,
                               postprocess);

  /*! [mesh_warping] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mesh smoothing.
 *
 * \param[in,out] mesh pointer to cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_smoothe(cs_mesh_t  *mesh)
{
  /*! [mesh_smoothing] */

  double feature_angle = 25; /* bounded between 0 and 90 degrees */
  int *vtx_is_fixed = NULL;

  BFT_MALLOC(vtx_is_fixed, mesh->n_vertices, int);

  /* Get fixed boundary vertices flag */

  cs_mesh_smoother_fix_by_feature(mesh,
                                  feature_angle,
                                  vtx_is_fixed);

  /* Call unwarping smoother */

  cs_mesh_smoother_unwarp(mesh, vtx_is_fixed);

  /* Free memory */

  BFT_FREE(vtx_is_fixed);

  /*! [mesh_smoothing] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tag bad cells within the mesh based on user-defined geometric criteria.
 *
 * \param[in,out] mesh pointer to cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to cs_mesh_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_bad_cells_tag(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities)
{
  /* Example: tag cells having a volume below 0.01 m^3 */
  /*          and post-process the tagged cells        */
  /*---------------------------------------------------*/

  /*! [mesh_tag_bad_cells] */

  const cs_lnum_t  n_cells         = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;
  unsigned  *bad_cell_flag         = mesh_quantities->bad_cell_flag;

  double *volume  = mesh_quantities->cell_vol;

  cs_lnum_t cell_id;
  cs_gnum_t n_cells_tot, iwarning, ibad;

  /*---------------------------------------*/
  /* User condition: check for cell volume */
  /*---------------------------------------*/

  cs_lnum_t  *bad_vol_cells = NULL;

  BFT_MALLOC(bad_vol_cells, n_cells_wghosts, cs_lnum_t);

  for (cell_id = 0; cell_id < n_cells_wghosts; cell_id++)
    bad_vol_cells[cell_id] = 0;

  /* Get the total number of cells within the domain */
  n_cells_tot = n_cells;
  cs_parall_counter(&n_cells_tot, 1);

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Compare the cell volume to the user condition */

    if (volume[cell_id] < 0.01) {
      /* Local array used to post-process results
         --> user tagged cells are set to 1 */
      bad_vol_cells[cell_id] = 1;

      /* Array used to store bad cells flag
         --> user tagged cells are flagged using the mask */
      bad_cell_flag[cell_id] = bad_cell_flag[cell_id] | CS_BAD_CELL_USER;
    }
  }

  ibad = 0;
  iwarning = 0;
  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    if (bad_cell_flag[cell_id] & CS_BAD_CELL_USER) {
      ibad++;
      iwarning++;
    }
  }

  /* Parallel processing */
  if (cs_glob_rank_id >= 0) {
    cs_parall_counter(&ibad, 1);
    cs_parall_counter(&iwarning, 1);
  }

  /* Display log output */
  bft_printf("\n  Criteria 6: User Specific Tag:\n");
  bft_printf("    Number of bad cells detected: %llu --> %3.0f %%\n",
             (unsigned long long)ibad,
             (double)ibad / (double)n_cells_tot * 100.0);

  if (iwarning > 0)
    bft_printf
      (" Warning:\n"
       " --------\n"
       "    Mesh quality issue based on user criteria has been detected\n\n"
       "    The mesh should be re-considered using the listed criteria.\n\n");

  /* Post processing is activated automatically in mesh quality mode
     for the CS_BAD_CELL_USER tag, as for all tags defined in
     cs_user_bad_cells.h.
     For a different tag, postprocessing could be done with a user
     postprocessing function, or calling directly cs_post_write_var(). */

  BFT_FREE(bad_vol_cells);

  /*! [mesh_tag_bad_cells] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
