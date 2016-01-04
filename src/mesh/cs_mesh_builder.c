/*============================================================================
 * \file Auxiliary structure used to read, write, and partition mesh data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_periodicity.h"
#include "fvm_selector.h"

#include "cs_base.h"
#include "cs_interface.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the temporary mesh used for building main mesh */

cs_mesh_builder_t  *cs_glob_mesh_builder = NULL;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an empty mesh builder structure.
 *
 * \return  pointer to a mesh builder structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_create(void)
{
  cs_mesh_builder_t  *mb = NULL;

  BFT_MALLOC(mb, 1, cs_mesh_builder_t);

  mb->n_g_faces = 0;
  mb->n_g_face_connect_size = 0;

  mb->n_perio = 0;

  mb->have_cell_rank = false;

  /* Temporary mesh data */

  mb->face_cells = NULL;
  mb->face_vertices_idx = NULL;
  mb->face_vertices = NULL;
  mb->cell_gc_id = NULL;
  mb->face_gc_id = NULL;
  mb->vertex_coords = NULL;

  /* Periodic features */

  mb->periodicity_num = NULL;
  mb->n_per_face_couples = NULL;
  mb->n_per_face_couples = NULL;
  mb->n_g_per_face_couples = NULL;
  mb->per_face_couples = NULL;

  /* Optional partitioning info */

  mb->cell_rank = NULL;

  /* Block ranges for parallel distribution */

  mb->min_rank_step = 1;
  memset(&(mb->cell_bi), 0, sizeof(cs_block_dist_info_t));
  memset(&(mb->face_bi), 0, sizeof(cs_block_dist_info_t));
  memset(&(mb->vertex_bi), 0, sizeof(cs_block_dist_info_t));
  mb->per_face_bi = NULL;

  return mb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_mesh_builder_t structure.
 *
 * \param[in, out]  mb  pointer to pointer of structure to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_builder_destroy(cs_mesh_builder_t  **mb)
{
  if (mb != NULL) {

    cs_mesh_builder_t  *_mb = *mb;

    /* Temporary mesh data */

    BFT_FREE(_mb->face_cells);
    BFT_FREE(_mb->face_vertices_idx);
    BFT_FREE(_mb->face_vertices);
    BFT_FREE(_mb->cell_gc_id);
    BFT_FREE(_mb->face_gc_id);
    BFT_FREE(_mb->vertex_coords);

    /* Periodic features */

    BFT_FREE(_mb->periodicity_num);
    BFT_FREE(_mb->n_per_face_couples);
    BFT_FREE(_mb->n_g_per_face_couples);
    if (_mb->per_face_couples != NULL) {
      for (int i = 0; i < _mb->n_perio; i++)
        BFT_FREE(_mb->per_face_couples[i]);
      BFT_FREE(_mb->per_face_couples);
    }

    /* Optional partitioning info */

    BFT_FREE(_mb->cell_rank);

    /* Block ranges for parallel distribution */

    BFT_FREE(_mb->per_face_bi);

    BFT_FREE(*mb);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define block distribution sizes for mesh builder.
 *
 * \param[in, out]  mb              pointer to mesh builder to update
 * \param[in]       rank_id         id of local rank
 * \param[in]       n_ranks         number of associated ranks
 * \param[in]       min_rank_step   minimum rank step between blocks
 * \param[in]       min_block_size  minimum number of entities per block
 * \param[in]       n_g_cells       global number of cells
 * \param[in]       n_g_faces       global number of faces
 * \param[in]       n_g_vertices    global number of vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_builder_define_block_dist(cs_mesh_builder_t  *mb,
                                  int                 rank_id,
                                  int                 n_ranks,
                                  int                 min_rank_step,
                                  int                 min_block_size,
                                  cs_gnum_t           n_g_cells,
                                  cs_gnum_t           n_g_faces,
                                  cs_gnum_t           n_g_vertices)
{
  mb->min_rank_step = min_rank_step;

  mb->cell_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            min_rank_step,
                                            min_block_size,
                                            n_g_cells);

  mb->face_bi = cs_block_dist_compute_sizes(rank_id,
                                            n_ranks,
                                            min_rank_step,
                                            min_block_size,
                                            n_g_faces);

  mb->vertex_bi = cs_block_dist_compute_sizes(rank_id,
                                              n_ranks,
                                              min_rank_step,
                                              min_block_size,
                                              n_g_vertices);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
