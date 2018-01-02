/*============================================================================
 * \file Auxiliary structure used to read, write, and partition mesh data.
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

#include "cs_mesh_builder.h"

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

/*=============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compare periodic couples in global numbering form (qsort function).
 *
 * parameters:
 *   x <-> pointer to first couple
 *   y <-> pointer to second couple
 *
 * returns:
 *   lexicographical
 *----------------------------------------------------------------------------*/

static int _compare_couples(const void *x, const void *y)
{
  int retval = 1;

  const cs_gnum_t *c0 = x;
  const cs_gnum_t *c1 = y;

  if (c0[0] < c1[0])
    retval = -1;

  else if (c0[0] == c1[0]) {
    if (c0[1] < c1[1])
      retval = -1;
    else if (c0[1] == c1[1])
      retval = 0;
  }

  return retval;
}

#endif /* defined(HAVE_MPI) */

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

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract periodic face connectivity information from faces interface
 *        for mesh builder when running in parallel mode.
 *
 * \param[in]       n_init_perio  number of initial periodicities
 * \param[in, out]  mb            pointer to mesh builder structure
 * \param[in]       periodicity   periodicity information
 * \param[in]       face_gnum     global face numbers, or NULL
 * \param[in]       face_ifs      parallel and periodic faces interfaces set
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_builder_extract_periodic_faces_g(int                        n_init_perio,
                                         cs_mesh_builder_t         *mb,
                                         fvm_periodicity_t         *periodicity,
                                         const cs_gnum_t           *face_gnum,
                                         const cs_interface_set_t  *face_ifs)
{
  int i, j;
  cs_lnum_t k, l;

  int perio_count = 0;
  cs_lnum_t  *send_index = NULL;
  cs_gnum_t  *recv_num = NULL;
  int  *tr_id = NULL;

  cs_datatype_t gnum_type = CS_GNUM_TYPE;

  const int n_perio = n_init_perio;
  const int n_interfaces = cs_interface_set_size(face_ifs);

  /* Free previous values if we are updating */

  if (mb->n_perio > 0 && mb->n_per_face_couples != NULL) {
    for (i = 0; i < n_perio; i++)
      BFT_FREE(mb->per_face_couples[i]);
    BFT_FREE(mb->n_per_face_couples);
    BFT_FREE(mb->per_face_couples);
  }

  mb->n_perio = n_perio;

  /* Allocate arrays in mesh builder (initializing per_face_idx) */

  assert(periodicity != NULL);
  assert(mb != NULL);
  assert(mb->n_g_per_face_couples == 0);

  BFT_MALLOC(mb->n_per_face_couples, n_perio, cs_lnum_t);
  BFT_MALLOC(mb->per_face_couples, n_perio, cs_gnum_t *);

  for (i = 0; i < n_perio; i++) {
    mb->n_per_face_couples[i] = 0;
    mb->per_face_couples[i] = NULL;
  }

  /* List direct and reverse transforms */

  BFT_MALLOC(tr_id, n_perio*2, int);

  for (i = 0; i < n_perio*2; i++) {
    int rev_id = fvm_periodicity_get_reverse_id(periodicity, i);
    if (i < rev_id) {
      int parent_ids[2];
      fvm_periodicity_get_parent_ids(periodicity, i, parent_ids);
      if (parent_ids[0] < 0 && parent_ids[1] < 0) {
        tr_id[perio_count*2] = i + 1;
        tr_id[perio_count*2 + 1] = rev_id + 1;
        perio_count++;
      }
    }
  }
  assert(perio_count == n_perio);

  for (i = 0; i < n_interfaces; i++) {
    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, i);
    const cs_lnum_t *tr_index = cs_interface_get_tr_index(face_if);
    for (j = 0; j < n_perio; j++) {
      const cs_lnum_t n_tr_faces = (  tr_index[tr_id[j*2] + 1]
                                    - tr_index[tr_id[j*2]]);
      mb->n_per_face_couples[j] += n_tr_faces;
    }
  }

  BFT_MALLOC(recv_num, cs_interface_set_n_elts(face_ifs), cs_gnum_t);

  cs_interface_set_copy_array(face_ifs,
                              gnum_type,
                              1,
                              true, /* src_on_parent */
                              face_gnum,
                              recv_num);

  /* Prepare send buffer (send reverse transformation values) */

  BFT_FREE(send_index);

  for (i = 0; i < n_perio; i++)
    BFT_MALLOC(mb->per_face_couples[i], mb->n_per_face_couples[i]*2, cs_gnum_t);

  /* Reset couples count */

  for (i = 0; i < n_perio; i++)
    mb->n_per_face_couples[i] = 0;

  /* Copy face couples to mesh builder */

  for (i = 0, j = 0, l = 0; i < n_interfaces; i++) {

    const cs_interface_t *face_if = cs_interface_set_get(face_ifs, i);
    const cs_lnum_t *tr_index = cs_interface_get_tr_index(face_if);
    const cs_lnum_t *elt_id = cs_interface_get_elt_ids(face_if);

    l += tr_index[1];

    for (j = 0; j < n_perio; j++) {

      /* Count couples in direct periodicity */

      cs_lnum_t nc = mb->n_per_face_couples[j]*2;
      const cs_lnum_t start_id = tr_index[tr_id[j*2]];
      const cs_lnum_t end_id = tr_index[tr_id[j*2] + 1];

      for (k = start_id; k < end_id; k++) {
        cs_lnum_t f_id = elt_id[k];
        mb->per_face_couples[j][nc++] = face_gnum[f_id];
        mb->per_face_couples[j][nc++] = recv_num[l++];
      }
      mb->n_per_face_couples[j] = nc/2;

      /* Ignore couples in reverse periodicity */

      l += tr_index[tr_id[j*2 + 1] + 1] - tr_index[tr_id[j*2 + 1]];

    }

  }

  BFT_FREE(recv_num);
  BFT_FREE(tr_id);

  /* Now sort couples in place for future use (more for consistency
     and ease of verification than absolutely necessary) */

  for (i = 0; i < n_perio; i++) {
    if (mb->n_per_face_couples[i] > 0)
      qsort(mb->per_face_couples[i],
            mb->n_per_face_couples[i],
            sizeof(cs_gnum_t) * 2,
            &_compare_couples);
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
