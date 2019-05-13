/*============================================================================
 * Filters for dynamic models.
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

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_les_filter.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute filters for dynamic models.
 *
 * This function deals with the standard or extended neighborhood.
 *
 * \param[in]   stride   stride of array to filter
 * \param[in]   val      array of values to filter
 * \param[out]  f_val    array of filtered values
 */
/*----------------------------------------------------------------------------*/

void
cs_les_filter(int        stride,
              cs_real_t  val[],
              cs_real_t  f_val[])
{
  cs_real_t *w1 = NULL, *w2 = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_lnum_t  _stride = stride;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_elts_l = n_cells * _stride;
  const cs_lnum_t  n_elts = n_cells_ext * _stride;
  const int n_i_groups = mesh->i_face_numbering->n_groups;
  const int n_i_threads = mesh->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = mesh->i_face_numbering->group_index;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  assert(cell_cells_idx != NULL);

  /* Allocate and initialize working buffers */

  BFT_MALLOC(w1, n_elts, cs_real_t);
  BFT_MALLOC(w2, n_elts, cs_real_t);

  /* Case for scalar variable */
  /*--------------------------*/

  if (stride == 1) {

    /* Synchronize valiable */

    if (mesh->halo != NULL)
      cs_halo_sync_var(mesh->halo, CS_HALO_EXTENDED, val);

    /* Define filtered valiable array */

#   pragma omp parallel for
    for (cs_lnum_t i = 0; i < n_cells; i++) {

      w1[i] = val[i] * cell_vol[i];
      w2[i] = cell_vol[i];

      /* Loop on connected cells (without cells sharing a face) */

      for (cs_lnum_t j = cell_cells_idx[i]; j < cell_cells_idx[i+1]; j++) {
        cs_lnum_t k = cell_cells_lst[j];
        w1[i] += val[k] * cell_vol[k];
        w2[i] += cell_vol[k];
      }

    } /* End of loop on cells */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t i = mesh->i_face_cells[face_id][0];
          cs_lnum_t j = mesh->i_face_cells[face_id][1];

          w1[i] += val[j] * cell_vol[j];
          w2[i] += cell_vol[j];
          w1[j] += val[i] * cell_vol[i];
          w2[j] += cell_vol[i];

        }

      }

    }

#   pragma omp parallel for
    for (cs_lnum_t i = 0; i < n_cells; i++)
      f_val[i] = w1[i]/w2[i];

    /* Synchronize valiable */

    if (mesh->halo != NULL)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_val);
  }

  /* Case for vector or tensor variable */
  /*------------------------------------*/

  else {

    /* Synchronize valiable */

    if (mesh->halo != NULL)
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED, val, stride);

    /* Define filtered valiable array */

#   pragma omp parallel for
    for (cs_lnum_t i = 0; i < n_cells; i++) {

      for (cs_lnum_t c_id = 0; c_id < _stride; c_id++) {
        const cs_lnum_t ic = i *_stride + c_id;
        w1[ic] = val[ic] * cell_vol[i];
        w2[ic] = cell_vol[i];
      }

      /* Loop on connected cells (without cells sharing a face) */

      for (cs_lnum_t j = cell_cells_idx[i]; j < cell_cells_idx[i+1]; j++) {
        cs_lnum_t k = cell_cells_lst[j];
        for (cs_lnum_t c_id = 0; c_id < _stride; c_id++) {
          const cs_lnum_t ic = i *_stride + c_id;
          const cs_lnum_t kc = k *_stride + c_id;
          w1[ic] += val[kc] * cell_vol[k];
          w2[ic] += cell_vol[k];
        }
      }

    } /* End of loop on cells */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {

        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t i = mesh->i_face_cells[face_id][0];
          cs_lnum_t j = mesh->i_face_cells[face_id][1];

          for (cs_lnum_t c_id = 0; c_id < _stride; c_id++) {
            const cs_lnum_t ic = i *_stride + c_id;
            const cs_lnum_t jc = j *_stride + c_id;
            w1[ic] += val[jc] * cell_vol[j];
            w2[ic] += cell_vol[j];
            w1[jc] += val[ic] * cell_vol[i];
            w2[jc] += cell_vol[i];
          }

        }

      }

    }

#   pragma omp parallel for
    for (cs_lnum_t i = 0; i < n_elts_l; i++)
      f_val[i] = w1[i]/w2[i];

    /* Synchronize valiable */

    if (mesh->halo != NULL)
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED, f_val, stride);
  }

  BFT_FREE(w2);
  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
