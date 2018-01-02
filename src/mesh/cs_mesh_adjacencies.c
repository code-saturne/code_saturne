/*============================================================================
 * Additional mesh adjacencies.
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
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_adjacencies.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_mesh_adjacencies.c
        Additional mesh adjacencies.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static cs_mesh_adjacencies_t  _cs_glob_mesh_adjacencies;

const cs_mesh_adjacencies_t  *cs_glob_mesh_adjacencies = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update cell -> cells connectivity
 *
 * parameters:
 *   ma <-> mesh adjacecies structure to update
 *----------------------------------------------------------------------------*/

static void
_update_cell_cells(cs_mesh_adjacencies_t  *ma)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_2_t *restrict face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_faces = m->n_i_faces;

  /* Allocate and map */

  BFT_REALLOC(ma->cell_cells_idx, n_cells + 1, cs_lnum_t);
  cs_lnum_t *c2c_idx = ma->cell_cells_idx;

  /* Count number of nonzero elements per row */

  cs_lnum_t  *count;

  BFT_MALLOC(count, n_cells, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_cells; i++)
    count[i] = 0;

  for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
    cs_lnum_t i = face_cells[face_id][0];
    cs_lnum_t j = face_cells[face_id][1];
    if (i < n_cells)
      count[i] += 1;
    if (j < n_cells)
      count[j] += 1;
  }

  c2c_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    c2c_idx[i+1] = c2c_idx[i] + count[i];
    count[i] = 0;
  }

  /* Build structure */

  BFT_REALLOC(ma->cell_cells, c2c_idx[n_cells], cs_lnum_t);

  cs_lnum_t *c2c = ma->cell_cells;

  for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
    cs_lnum_t i = face_cells[face_id][0];
    cs_lnum_t j = face_cells[face_id][1];
    if (i < n_cells) {
      c2c[c2c_idx[i] + count[i]] = j;
      count[i] += 1;
    }
    if (j < n_cells) {
      c2c[c2c_idx[j] + count[j]] = i;
      count[j] += 1;
    }
  }

  BFT_FREE(count);

  /* Sort line elements by column id (for better access patterns) */

  ma->single_faces_to_cells = cs_sort_indexed(n_cells, c2c_idx, c2c);

  /* Compact elements if necessary */

  if (ma->single_faces_to_cells == false) {

    cs_lnum_t *tmp_c2c_idx = NULL;

    BFT_MALLOC(tmp_c2c_idx, n_cells+1, cs_lnum_t);
    memcpy(tmp_c2c_idx, c2c_idx, (n_cells+1)*sizeof(cs_lnum_t));

    cs_lnum_t k = 0;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      cs_lnum_t js = tmp_c2c_idx[i];
      cs_lnum_t je = tmp_c2c_idx[i+1];
      cs_lnum_t c2c_prev = -1;
      c2c_idx[i] = k;
      for (cs_lnum_t j = js; j < je; j++) {
        if (c2c_prev != c2c[j]) {
          c2c[k++] = c2c[j];
          c2c_prev = c2c[j];
        }
      }
    }
    c2c_idx[n_cells] = k;

    assert(c2c_idx[n_cells] < tmp_c2c_idx[n_cells]);

    BFT_FREE(tmp_c2c_idx);
    BFT_REALLOC(c2c, c2c_idx[n_cells], cs_lnum_t);

    ma->cell_cells = c2c;

  }
}

/*----------------------------------------------------------------------------
 * Update cells -> boundary faces connectivity
 *
 * parameters:
 *   ma <-> mesh adjacecies structure to update
 *----------------------------------------------------------------------------*/

static void
_update_cell_b_faces(cs_mesh_adjacencies_t  *ma)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* (re)build cell -> boundary faces index */

  BFT_REALLOC(ma->cell_b_faces_idx, n_cells + 1, cs_lnum_t);
  cs_lnum_t *c2b_idx = ma->cell_b_faces_idx;

  cs_lnum_t *c2b_count;
  BFT_MALLOC(c2b_count, n_cells, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_cells; i++)
    c2b_count[i] = 0;

  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    c2b_count[b_face_cells[i]] += 1;

  c2b_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    c2b_idx[i+1] = c2b_idx[i] + c2b_count[i];
    c2b_count[i] = 0;
  }

  /* Rebuild values */

  BFT_REALLOC(ma->cell_b_faces, c2b_idx[n_cells], cs_lnum_t);
  cs_lnum_t *c2b = ma->cell_b_faces;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t c_id = b_face_cells[i];
    c2b[c2b_idx[c_id] + c2b_count[c_id]] = i;
    c2b_count[c_id] += 1;
  }

  BFT_FREE(c2b_count);

  /* Sort array */

  cs_sort_indexed(n_cells, c2b_idx, c2b);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mesh adjacencies helper API.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_initialize(void)
{
  cs_mesh_adjacencies_t *ma = &_cs_glob_mesh_adjacencies;

  ma->cell_cells_idx = NULL;
  ma->cell_cells = NULL;

  ma->cell_cells_e_idx = NULL;
  ma->cell_cells_e = NULL;

  ma->cell_b_faces_idx = NULL;
  ma->cell_b_faces = NULL;

  cs_glob_mesh_adjacencies = ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize mesh adjacencies helper API.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_finalize(void)
{
  cs_mesh_adjacencies_t *ma = &_cs_glob_mesh_adjacencies;

  BFT_FREE(ma->cell_cells_idx);
  BFT_FREE(ma->cell_cells);

  BFT_FREE(ma->cell_b_faces_idx);
  BFT_FREE(ma->cell_b_faces);

  cs_glob_mesh_adjacencies = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update mesh adjacencies helper API relative to mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_mesh(void)
{
  cs_mesh_adjacencies_t *ma = &_cs_glob_mesh_adjacencies;

  /* (re)build cell -> cell connectivities */

  _update_cell_cells(ma);

  /* Map shared connectivities */

  cs_mesh_adjacencies_update_cell_cells_e();

  /* (re)build cell -> boundary face connectivities */

  _update_cell_b_faces(ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update extended cell -> cell connectivites in
 *         mesh adjacencies helper API relative to mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_adjacencies_update_cell_cells_e(void)
{
  const cs_mesh_t *m = cs_glob_mesh;

  cs_mesh_adjacencies_t *ma = &_cs_glob_mesh_adjacencies;

  ma->cell_cells_e_idx = m->cell_cells_idx;
  ma->cell_cells_e = m->cell_cells_lst;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
