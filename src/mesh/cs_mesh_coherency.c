/*============================================================================
 * Functions checking the coherency of the mesh
 *============================================================================*/

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

#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_coherency.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 *  Global static variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check the coherency of the internal face -> cells connectivity.
 *----------------------------------------------------------------------------*/

static void
_check_i_face_cells(void)
{
  cs_lnum_t  i;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  bft_printf(_("    Checking the face -> cells connectivity coherency\n"));

  for (i = 0; i < mesh->n_i_faces; i++) {

    if (mesh->i_face_cells[i][0] == -1 || mesh->i_face_cells[i][1] == -1)
      bft_error(__FILE__, __LINE__, 0,
                _("Internal face -> cells connectivity value not initialized\n"
                  "for face: %d (cell_num1 = %d and cell_num2 = %d)\n"),
                i+1, mesh->i_face_cells[i][0], mesh->i_face_cells[i][1]);

  }

}

/*----------------------------------------------------------------------------
 * Check that bounding boxes of 2 related cells intersect.
 *
 * parameters:
 *   halo_name <-- name of halo type
 *   cell_id1  <-- id of first cell
 *   cell_id2  <-- id of second cell
 *   emin      <-- minimum coordinates of bounding boxes for cells
 *   emax      <-- maximum coordinates of bounding boxes for cells
 *----------------------------------------------------------------------------*/

static void
_check_bounding_boxes(const char        *halo_type,
                      cs_lnum_t          cell_id1,
                      cs_lnum_t          cell_id2,
                      const cs_real_3_t  emin[],
                      const cs_real_3_t  emax[])
{
  int i;

  for (i = 0; i < 3; i++) {

    cs_real_t  delta1 = (emax[cell_id1][i] - emin[cell_id1][i])  * 0.5025;
    cs_real_t  delta2 = (emax[cell_id2][i] - emin[cell_id2][i])  * 0.5025;
    cs_real_t  mean1 = (emax[cell_id1][i] + emin[cell_id1][i])  * 0.5;
    cs_real_t  mean2 = (emax[cell_id2][i] + emin[cell_id2][i])  * 0.5;
    cs_real_t  bmin1 = mean1 - delta1;
    cs_real_t  bmax1 = mean1 + delta1;
    cs_real_t  bmin2 = mean2 - delta2;
    cs_real_t  bmax2 = mean2 + delta2;

    if (! (   (mean2 >= mean1 && bmin2 < bmax1)
           || (mean2 <  mean1 && bmin1 < bmax2))) {

      bft_error(__FILE__, __LINE__, 0,
                _("\nCoherency error in %s halo\n"
                  "between cell %ld with:\n"
                  "  bounding box min:  [%12.6g %12.6g %12.6g]\n"
                  "               max:  [%12.6g %12.6g %12.6g]\n"
                  "and     cell %ld with:\n"
                  "  bounding box min:  [%12.6g %12.6g %12.6g]\n"
                  "               max:  [%12.6g %12.6g %12.6g]"),
                halo_type,
                (long)cell_id1+1,
                emin[cell_id1][0], emin[cell_id1][1], emin[cell_id1][2],
                emax[cell_id1][0], emax[cell_id1][1], emax[cell_id1][2],
                (long)cell_id2+1,
                emin[cell_id2][0], emin[cell_id2][1], emin[cell_id2][2],
                emax[cell_id2][0], emax[cell_id2][1], emax[cell_id2][2]);

    }

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check the coherency of the global mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_coherency_check(void)
{
  cs_lnum_t  i, j, k, cell_id, face_id;

  cs_real_3_t  *emin = NULL, *emax = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const cs_lnum_t  *b_face_cells = mesh->b_face_cells;
  const cs_real_t  *vtx_coord = mesh->vtx_coord;

  bft_printf(_("\n Checking the mesh structure coherency:\n"));

  /* Check internal face -> cells connectivity coherency */

  _check_i_face_cells();

  /* allocate and initialize buffers */

  BFT_MALLOC(emin, n_cells_with_ghosts, cs_real_3_t);
  BFT_MALLOC(emax, n_cells_with_ghosts, cs_real_3_t);

  /* Loop on coordinates */

  bft_printf(_("    Coherency criteria definition\n"));

  for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    for (i = 0; i < 3; i++) {
      emin[cell_id][i] =  DBL_MAX;
      emax[cell_id][i] = -DBL_MAX;
    }
  }

  /* Loop on internal faces */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    const cs_lnum_t *face_vtx_idx = mesh->i_face_vtx_idx;
    const cs_lnum_t *face_vtx_lst = mesh->i_face_vtx_lst;

    cs_real_t _min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    cs_real_t _max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

    for (i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {

      cs_lnum_t vtx_id = face_vtx_lst[i];
      const cs_real_t *coord = vtx_coord + (3*vtx_id);

      for (j = 0; j < 3; j++) {
        _min[j] = CS_MIN(_min[j], coord[j]);
        _max[j] = CS_MAX(_max[j], coord[j]);
      }

    }

    for (j = 0; j < 2; j++) {

      cell_id = i_face_cells[face_id][j];

      for (k = 0; k < 3; k++) {
        emin[cell_id][k] = CS_MIN(emin[cell_id][k], _min[k]);
        emax[cell_id][k] = CS_MAX(emax[cell_id][k], _max[k]);
      }

    }

  } /* End of loop on internal faces */

  /* Loop on border faces */

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {

    const cs_lnum_t *face_vtx_idx = mesh->b_face_vtx_idx;
    const cs_lnum_t *face_vtx_lst = mesh->b_face_vtx_lst;

    cs_real_t _min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    cs_real_t _max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

    for (i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {

      cs_lnum_t vtx_id = face_vtx_lst[i];
      const cs_real_t *coord = vtx_coord + (3*vtx_id);

      for (j = 0; j < 3; j++) {
        _min[j] = CS_MIN(_min[j], coord[j]);
        _max[j] = CS_MAX(_max[j], coord[j]);
      }

      cell_id = b_face_cells[face_id];

      for (j = 0; j < 3; j++) {
        emin[cell_id][j] = CS_MIN(emin[cell_id][j], _min[j]);
        emax[cell_id][j] = CS_MAX(emax[cell_id][j], _max[j]);
      }

    }

  } /* End of loop on border faces */

  /* Synchronize variable */

  if (mesh->halo != NULL) {
    cs_halo_sync_var_strided(mesh->halo, mesh->halo_type, (cs_real_t *)emin, 3);
    cs_halo_sync_var_strided(mesh->halo, mesh->halo_type, (cs_real_t *)emax, 3);
  }

  /* Synchronization for periodicity */

  if (mesh->n_init_perio > 0) {
    cs_halo_perio_sync_coords(mesh->halo, mesh->halo_type, (cs_real_t *)emin);
    cs_halo_perio_sync_coords(mesh->halo, mesh->halo_type, (cs_real_t *)emax);
  }

  /* For rotational periodicity, things are more complex, as
     rotation of a bounding box is implied */

  if (mesh->have_rotation_perio) {

    int corner;

    cs_real_3_t  *smin = NULL, *smax = NULL;
    cs_real_3_t  *c_coords = NULL;

    /* allocate and initialize buffers */

    BFT_MALLOC(smin, n_cells_with_ghosts, cs_real_3_t);
    BFT_MALLOC(smax, n_cells_with_ghosts, cs_real_3_t);
    BFT_MALLOC(c_coords, n_cells_with_ghosts, cs_real_3_t);

    memcpy(smin, emin, n_cells_with_ghosts*sizeof(cs_real_3_t));
    memcpy(smax, emax, n_cells_with_ghosts*sizeof(cs_real_3_t));

    /* Reset min/max in halo */

    for (cell_id = mesh->n_cells; cell_id < n_cells_with_ghosts; cell_id++) {
      for (i = 0; i < 3; i++) {
        emin[cell_id][i] =  DBL_MAX;
        emax[cell_id][i] = -DBL_MAX;
      }
    }

    for (corner = 0; corner < 8; corner++) {

      cs_real_3_t *px = (((corner+1) % 4) < 2) ? smin : smax;
      cs_real_3_t *py = ((corner % 4) < 2) ? smin : smax;
      cs_real_3_t *pz = (corner < 4) ? smin : smax;

      for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
        c_coords[cell_id][0] = px[cell_id][0];
        c_coords[cell_id][1] = py[cell_id][1];
        c_coords[cell_id][2] = pz[cell_id][2];
      }

      cs_halo_sync_var_strided
        (mesh->halo, mesh->halo_type, (cs_real_t *)c_coords, 3);
      cs_halo_perio_sync_coords
        (mesh->halo, mesh->halo_type, (cs_real_t *)c_coords);

      for (cell_id = mesh->n_cells; cell_id < n_cells_with_ghosts; cell_id++) {
        for (j = 0; j < 3; j++) {
          emin[cell_id][j] = CS_MIN(emin[cell_id][j], c_coords[cell_id][j]);
          emax[cell_id][j] = CS_MAX(emax[cell_id][j], c_coords[cell_id][j]);
        }
      }

    }

    BFT_FREE(smin);
    BFT_FREE(smax);
    BFT_FREE(c_coords);

  } /* end of additional updates for rotational periodicity */

  bft_printf(_("    Coherency verification on coordinates\n"));

  /* Test coherency on the standard neighborhood */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cs_lnum_t  cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t  cell_id2 = i_face_cells[face_id][1];

    _check_bounding_boxes(_("standard"),
                          cell_id1,
                          cell_id2,
                          (const cs_real_3_t *)emin,
                          (const cs_real_3_t *)emax);

  } /* End of loop on internal faces */

  if (mesh->cell_cells_idx != NULL) {

    cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      for (i = cell_cells_idx[cell_id];
           i < cell_cells_idx[cell_id+1];
           i++) {

        cs_lnum_t  cell_id2 = mesh->cell_cells_lst[i];

        _check_bounding_boxes(_("extended"),
                              cell_id,
                              cell_id2,
                              (const cs_real_3_t *)emin,
                              (const cs_real_3_t *)emax);

      }

    } /* End of loop on cells */

  } /* End of treatment of the exetended neighborhood */

  /* Free memory */

  BFT_FREE(emin);
  BFT_FREE(emax);

  bft_printf(_(" End of coherency check of the mesh structure.\n"));

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
