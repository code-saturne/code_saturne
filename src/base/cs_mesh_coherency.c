/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Functions checking the coherency of the mesh
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_perio.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_coherency.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure and type definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_MESH_COHERENCY_TOLERANCE  0.05

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
_check_ifacel(void)
{
  cs_int_t  i;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  bft_printf(_("    Checking the face -> cells connectivity coherency\n"));

  for (i = 0; i < mesh->n_i_faces; i++) {

    if (mesh->i_face_cells[2*i] == 0 || mesh->i_face_cells[2*i+1] == 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Internal face -> cells connectivity value not initialized\n"
                  "for face: %d (cell_num1 = %d and cell_num2 = %d)\n"),
                i+1, mesh->i_face_cells[2*i], mesh->i_face_cells[2*i+1]);

  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check the coherency of the global mesh structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_coherency_check(void)
{
  cs_int_t  i, j, coord_id, cell_id, face_id, vtx_id;
  cs_real_t  test, delta_mean, delta_neighbor;
  cs_real_t  _min, _max, coord;

  cs_real_t delta_mean_mult = 1.0;

  cs_real_t  *minmax_buffer = NULL, *compute_buffer = NULL;
  cs_real_t  *min = NULL, *max = NULL, *mean = NULL, *delta = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;
  const cs_int_t  n_cells = mesh->n_cells;
  const cs_int_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  const cs_int_t  *ifacel = mesh->i_face_cells;
  const cs_int_t  *bfacel = mesh->b_face_cells;
  const cs_real_t  *vtx_coord = mesh->vtx_coord;

  bft_printf(_("\n Checking the mesh structure coherency:\n"));

  /* Check internal face -> cells connectivity coherency */

  _check_ifacel();

  /* allocate and initialize buffers */

  BFT_MALLOC(minmax_buffer, 2*n_cells, cs_real_t);
  BFT_MALLOC(compute_buffer, 6*n_cells_with_ghosts, cs_real_t);

  min = minmax_buffer;
  max = minmax_buffer + n_cells;
  mean = compute_buffer;
  delta = compute_buffer + 3*n_cells_with_ghosts;

  /* Loop on coordinates */

  bft_printf(_("    Coherency criteria definition\n"));

  for (coord_id = 0; coord_id < 3; coord_id++) {

    for (i = 0; i < n_cells; i++) {
      min[i] = DBL_MAX;
      max[i] =-DBL_MAX;
    }

    for (i = 0; i < n_cells_with_ghosts; i++) {
      mean[3*i + coord_id] = DBL_MIN;
      delta[3*i + coord_id] = DBL_MIN;
    }

    /* Loop on internal faces */

    for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

      const cs_int_t *face_vtx_idx = mesh->i_face_vtx_idx;
      const cs_int_t *face_vtx_lst = mesh->i_face_vtx_lst;

      _min = DBL_MAX;
      _max =-DBL_MAX;

      for (i = face_vtx_idx[face_id]-1; i < face_vtx_idx[face_id+1]-1; i++) {

        vtx_id = face_vtx_lst[i]-1;
        coord = vtx_coord[3*vtx_id + coord_id];
        _min = CS_MIN(_min, coord);
        _max = CS_MAX(_max, coord);

      }

      for (j = 0; j < 2; j++) {

        cell_id = ifacel[2*face_id + j] - 1;

        if (cell_id < n_cells) {
          max[cell_id] = CS_MAX(max[cell_id], _max);
          min[cell_id] = CS_MIN(min[cell_id], _min);
        }

      }

    } /* End of loop on internal faces */

    /* Loop on border faces */

    for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {

      const cs_int_t *face_vtx_idx = mesh->b_face_vtx_idx;
      const cs_int_t *face_vtx_lst = mesh->b_face_vtx_lst;

      _min = DBL_MAX;
      _max =-DBL_MAX;

      for (i = face_vtx_idx[face_id]-1; i < face_vtx_idx[face_id+1]-1; i++) {

        vtx_id = face_vtx_lst[i]-1;
        coord = vtx_coord[3*vtx_id + coord_id];
        _min = CS_MIN(_min, coord);
        _max = CS_MAX(_max, coord);

      }

      cell_id = bfacel[face_id] - 1;

      assert(cell_id < n_cells);
      max[cell_id] = CS_MAX(max[cell_id], _max);
      min[cell_id] = CS_MIN(min[cell_id], _min);

    } /* End of loop on border faces */

    /* Compute mean and delta for this coordinate */

    for (cell_id = 0; cell_id < n_cells; cell_id++) {

      mean[3*cell_id + coord_id] = (max[cell_id] + min[cell_id])/2;
      delta[3*cell_id + coord_id] = max[cell_id] - min[cell_id];

      assert(delta[3*cell_id + coord_id] > 0);

    }

  } /* End of loop on coordinates */

  /* Synchronize variable */

  if (mesh->halo != NULL) {

    cs_halo_sync_var_strided(mesh->halo, mesh->halo_type, mean, 3);
    cs_halo_sync_var_strided(mesh->halo, mesh->halo_type, delta, 3);

  }

  /* Synchronization for periodicity */

  if (mesh->n_init_perio > 0) {

    cs_real_t *delta_buffer;

    BFT_MALLOC(delta_buffer, 3*n_cells_with_ghosts, cs_real_t);

    /* De-interlace delta arrays for periodicity exchange;
       Also add factor of 1.8 (> sqrt(3)) as the longest diagonal
       of a box of side 1 is sqrt(3), and may find itself aligned
       with axes after rotation. */

    if (mesh->have_rotation_perio == 1) {

      delta_mean_mult = 1.0/1.8;

      for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
        cs_real_t delta_max = delta[3*cell_id];
        if (delta[3*cell_id + 1] > delta_max)
          delta_max = delta[3*cell_id + 1];
        if (delta[3*cell_id + 2] > delta_max)
          delta_max = delta[3*cell_id + 2];
        delta_max *= 1.8;
        delta_buffer[                        cell_id] = delta_max;
        delta_buffer[  n_cells_with_ghosts + cell_id] = delta_max;
        delta_buffer[2*n_cells_with_ghosts + cell_id] = delta_max;
      }

    }
    else {

      for (coord_id = 0; coord_id < 3; coord_id++) {
        for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++)
          delta_buffer[coord_id*n_cells_with_ghosts + cell_id]
            = delta[3*cell_id + coord_id];
      }

    }

    cs_perio_sync_var_vect(mesh->halo,
                           mesh->halo_type,
                           CS_PERIO_ROTA_COPY,
                           delta_buffer,
                           delta_buffer +   n_cells_with_ghosts,
                           delta_buffer + 2*n_cells_with_ghosts);

    for (coord_id = 0; coord_id < 3; coord_id++) {
      for (cell_id = n_cells; cell_id < n_cells_with_ghosts; cell_id++)
          delta[3*cell_id + coord_id] =
            delta_buffer[coord_id*n_cells_with_ghosts + cell_id];
    }

    BFT_FREE(delta_buffer);

    cs_perio_sync_coords(mesh->halo, mesh->halo_type, mean);

  }

  for (coord_id = 0; coord_id < 3; coord_id++) {

    bft_printf(_("    Coherency verification on coordinates %d\n"),
               coord_id+1);

    /* Test coherency on the standard neighborhood */

    for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

      cs_int_t  cell_id1 = ifacel[2*face_id] - 1;
      cs_int_t  cell_id2 = ifacel[2*face_id + 1] - 1;
      cs_real_t  delta1 = CS_ABS(delta[3*cell_id1 + coord_id]);
      cs_real_t  delta2 = CS_ABS(delta[3*cell_id2 + coord_id]);
      cs_real_t  mean1 = mean[3*cell_id1 + coord_id];
      cs_real_t  mean2 = mean[3*cell_id2 + coord_id];

      delta_mean = CS_ABS(mean2 - mean1)*delta_mean_mult;
      delta_neighbor = (delta1 + delta2)/2.;

      test = (1 + CS_MESH_COHERENCY_TOLERANCE)*delta_neighbor - delta_mean;

      if (test < 0) {

        cs_real_t  *cell_cen = mesh_quantities->cell_cen;

        bft_printf(_("\nInfo on cell 1: %d\n"
                     " cell center: %12.3g %12.3g %12.3g\n"
                     " delta:       %12.3g\n"
                     " box center:  %12.3g\n"),
                   cell_id1+1, cell_cen[3*cell_id1], cell_cen[3*cell_id1+1],
                   cell_cen[3*cell_id1+2], delta1, mean1);

        bft_printf(_("\nInfo on cell 2: %d\n"
                     " cell center: %12.3g %12.3g %12.3g\n"
                     " delta:       %12.3g\n"
                     " box center:  %12.3g\n"),
                   cell_id2+1, cell_cen[3*cell_id2], cell_cen[3*cell_id2+1],
                   cell_cen[3*cell_id2+2], delta2, mean2);
        bft_printf_flush();

        bft_error(__FILE__, __LINE__, 0,
                  _("\nCoherency error in standard halo\n"
                    "between cells %d and %d: test = %g\n"
                    "(delta = %g, delta_mean = %g)\n"),
                  cell_id1+1, cell_id2+1, test, delta_neighbor, delta_mean);

      }

    } /* End of loop on internal faces */

    if (mesh->cell_cells_idx != NULL) {

      cs_int_t  *cell_cells_idx = mesh->cell_cells_idx;

      for (cell_id = 0; cell_id < n_cells; cell_id++) {

        for (i = cell_cells_idx[cell_id]-1;
             i < cell_cells_idx[cell_id+1]-1; i++) {

          cs_int_t  cell_id2 = mesh->cell_cells_lst[i] - 1;
          cs_real_t  delta1 = CS_ABS(delta[3*cell_id + coord_id]);
          cs_real_t  delta2 = CS_ABS(delta[3*cell_id2 + coord_id]);
          cs_real_t  mean1 = mean[3*cell_id + coord_id];
          cs_real_t  mean2 = mean[3*cell_id2 + coord_id];

          delta_mean = CS_ABS(mean2 - mean1)*delta_mean_mult;
          delta_neighbor = (delta1 + delta2)/2.;

          test = (1 + CS_MESH_COHERENCY_TOLERANCE)*delta_neighbor - delta_mean;

          if (test < 0) {

            cs_real_t  *cell_cen = mesh_quantities->cell_cen;

            bft_printf(_("\nInfo on cell 1: %d\n"
                         " cell center: %12.3g %12.3g %12.3g\n"
                         " delta:       %12.3g\n"
                         " box center:  %12.3g\n"),
                       cell_id+1, cell_cen[3*cell_id], cell_cen[3*cell_id+1],
                       cell_cen[3*cell_id+2], delta1, mean1);

            bft_printf(_("\nInfo on cell 2: %d\n"
                         " cell center: %12.3g %12.3g %12.3g\n"
                         " delta:       %12.3g\n"
                         " box center:  %12.3g\n"),
                       cell_id2+1, cell_cen[3*cell_id2], cell_cen[3*cell_id2+1],
                       cell_cen[3*cell_id2+2], delta2, mean2);
            bft_printf_flush();

            bft_error(__FILE__, __LINE__, 0,
                      _("\nCoherency error in extended halo\n"
                        "between cells %d and %d: test = %g\n"
                        "(delta = %g, delta_mean = %g)\n"),
                      cell_id+1, cell_id2+1, test, delta_neighbor, delta_mean);

          }

        }

      } /* End of loop on cells */

    } /* End of treatment of the exetended neighborhood */

  } /* End of loop on coordinates */

  /* Free memory */

  BFT_FREE(compute_buffer);
  BFT_FREE(minmax_buffer);

  bft_printf(_(" End of coherency check of the mesh structure.\n"));

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
