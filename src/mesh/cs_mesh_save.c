/*============================================================================
 * Save mesh Preprocessor data
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_block_to_part.h"
#include "fvm_part_to_block.h"
#include "fvm_io_num.h"
#include "fvm_parall.h"
#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_order.h"
#include "cs_io.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_save.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write mesh periodicity metadata.
 *
 * parameters:
 *   mesh      <-- pointer to mesh structure
 *   perio_num <-- periodicity num
 *   pp_out    <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_perio_metadata(const cs_mesh_t  *mesh,
                           int               perio_num,
                           cs_io_t          *pp_out)
{
  char section_name[32];

  const cs_datatype_t lnum_type
    = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;

  double  matrix[3][4];
  cs_lnum_t perio_type = 0;

  const int tr_id = (perio_num - 1)*2;
  const fvm_periodicity_t *perio = mesh->periodicity;

  assert(perio != NULL);

  /* Get periodicity type and matrix */

  perio_type = fvm_periodicity_get_type(perio, tr_id);

  fvm_periodicity_get_matrix(perio, tr_id, matrix);

  /* Write global data */

  sprintf(section_name, "periodicity_type_%02d", perio_num);

  cs_io_write_global(section_name, 1, 0, 0, 1, lnum_type,
                     &perio_type,  pp_out);

  sprintf(section_name, "periodicity_matrix_%02d", perio_num);

  cs_io_write_global(section_name, 12, 0, 0, 1, CS_DOUBLE,
                     matrix,  pp_out);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write mesh data in parallel.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   face_bi  <-- face part to block info structure
 *   d        <-- face part to block distribution helper
 *   pp_out   <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_g(const cs_mesh_t           *mesh,
                       fvm_part_to_block_info_t   face_bi,
                       fvm_part_to_block_t       *d,
                       cs_io_t                   *pp_out)
{
  cs_lnum_t i, j, k;
  cs_gnum_t block_size, g_vtx_connect_size, n_block_faces, n_face_vertices;
  cs_gnum_t idx_range[4];

  cs_lnum_t *face_vtx_idx = NULL, *_face_vtx_idx = NULL;
  cs_gnum_t *face_vtx_idx_g = NULL;
  cs_gnum_t *face_vtx_g = NULL, *_face_vtx_g = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  BFT_MALLOC(face_vtx_idx, n_faces + 1, cs_lnum_t);
  BFT_MALLOC(_face_vtx_idx,
             (face_bi.gnum_range[1] - face_bi.gnum_range[0]) + 1,
             cs_lnum_t);

  face_vtx_idx[0] = 0;
  for (i = 0; i < n_i_faces; i++) {
    n_face_vertices = mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i];
    face_vtx_idx[i+1] = face_vtx_idx[i] + n_face_vertices;
  }
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
    n_face_vertices = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
    face_vtx_idx[j+1] = face_vtx_idx[j] + n_face_vertices;
  }

  fvm_part_to_block_copy_index(d, face_vtx_idx, _face_vtx_idx);

  /* Copy index from block to global values */

  n_block_faces = (face_bi.gnum_range[1] - face_bi.gnum_range[0]);
  BFT_MALLOC(face_vtx_idx_g, n_block_faces+ 1, cs_gnum_t);

  block_size = _face_vtx_idx[n_block_faces];
  MPI_Scan(&block_size, face_vtx_idx_g, 1, CS_MPI_GNUM, MPI_SUM,
           cs_glob_mpi_comm);
  face_vtx_idx_g[0] -= block_size;
  face_vtx_idx_g[0] += 1;

  for (i = 0; i < (cs_lnum_t)n_block_faces; i++) {
    n_face_vertices = _face_vtx_idx[i+1] - _face_vtx_idx[i];
    face_vtx_idx_g[i+1] = face_vtx_idx_g[i] + n_face_vertices;
  }

  idx_range[0] = face_bi.gnum_range[0];
  idx_range[1] = face_bi.gnum_range[1];
  if (face_bi.gnum_range[0] >= n_g_faces) {
    idx_range[0] += 1;
    idx_range[1] += 1;
  }
  else if (face_bi.gnum_range[1] >= n_g_faces + 1)
    idx_range[1] += 1;

  /* Set idx_range and compute global size for next write
     with indexed values */

  idx_range[2] = face_vtx_idx_g[0];
  idx_range[3] = face_vtx_idx_g[0] + block_size;

  MPI_Allreduce(face_vtx_idx_g + n_block_faces, &g_vtx_connect_size, 1,
                CS_MPI_GNUM, MPI_MAX, cs_glob_mpi_comm);
  g_vtx_connect_size -= 1;

  /* Now write buffer */

  cs_io_write_block_buffer("face_vertices_index",
                           n_g_faces + 1,
                           idx_range[0],
                           idx_range[1],
                           2, /* location_id, */
                           1, /* index id */
                           1, /* n_location_vals */
                           gnum_type,
                           face_vtx_idx_g,
                           pp_out);

  BFT_FREE(face_vtx_idx_g);

  /* Build connectivity */

  BFT_MALLOC(face_vtx_g,
             mesh->i_face_vtx_connect_size + mesh->b_face_vtx_connect_size,
             cs_gnum_t);

  k = 0;
  for (i = 0; i < n_i_faces; i++) {
    for (j = mesh->i_face_vtx_idx[i]; j < mesh->i_face_vtx_idx[i+1]; j++)
      face_vtx_g[k++]
        = mesh->global_vtx_num[mesh->i_face_vtx_lst[j - 1] - 1];
  }
  for (i = 0; i < n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1]; j++)
      face_vtx_g[k++]
        = mesh->global_vtx_num[mesh->b_face_vtx_lst[j - 1] - 1];
  }

  BFT_MALLOC(_face_vtx_g, block_size, cs_gnum_t);

  fvm_part_to_block_copy_indexed(d,
                                 gnum_type,
                                 face_vtx_idx,
                                 face_vtx_g,
                                 _face_vtx_idx,
                                 _face_vtx_g);

  BFT_FREE(face_vtx_g);
  BFT_FREE(_face_vtx_idx);
  BFT_FREE(face_vtx_idx);

  cs_io_write_block_buffer("face_vertices",
                           g_vtx_connect_size,
                           idx_range[2],
                           idx_range[3],
                           2, /* location_id, */
                           1, /* index id */
                           1, /* n_location_vals */
                           gnum_type,
                           _face_vtx_g,
                           pp_out);

  BFT_FREE(_face_vtx_g);
}

/*----------------------------------------------------------------------------
 * Write mesh data in parallel.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   pp_out   <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_data_g(const cs_mesh_t  *mesh,
                   cs_io_t          *pp_out)
{
  cs_lnum_t i, j;

  fvm_part_to_block_info_t cell_bi;
  fvm_part_to_block_info_t face_bi;
  fvm_part_to_block_info_t vtx_bi;
  fvm_part_to_block_t *d = NULL;

  cs_lnum_t *_cell_gc_id = NULL, *_face_gc_id = NULL;
  cs_lnum_t *face_gc_id = NULL;
  cs_gnum_t *cell_gnum = NULL, *face_gnum = NULL;
  cs_gnum_t *face_cell_g = NULL, *_face_cell_g = NULL;
  cs_real_t *_vtx_coords = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t real_type
    = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t lnum_type
    = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;

  /* Distribute cell group class info to blocks (write later) */
  /*----------------------------------------------------------*/

  cell_bi = fvm_part_to_block_compute_sizes(cs_glob_rank_id,
                                            cs_glob_n_ranks,
                                            0,
                                            0,
                                            mesh->n_g_cells);

  BFT_MALLOC(_cell_gc_id,
             (cell_bi.gnum_range[1] - cell_bi.gnum_range[0]),
             cs_lnum_t);

  d = fvm_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                       cell_bi,
                                       mesh->n_cells,
                                       mesh->global_cell_num);

  fvm_part_to_block_copy_array(d,
                               lnum_type,
                               1,
                               mesh->cell_family,
                               _cell_gc_id);

  fvm_part_to_block_destroy(&d);

  /* Build global face part to block distribution structures */
  /*---------------------------------------------------------*/

  BFT_MALLOC(face_gnum, n_faces, cs_gnum_t);

  for (i = 0; i < n_i_faces; i++)
    face_gnum[i] = mesh->global_i_face_num[i];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    face_gnum[j] = mesh->global_b_face_num[i] + mesh->n_g_i_faces;

  face_bi = fvm_part_to_block_compute_sizes(cs_glob_rank_id,
                                            cs_glob_n_ranks,
                                            0,
                                            0,
                                            n_g_faces);

  d = fvm_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                       face_bi,
                                       n_faces,
                                       face_gnum);

  /* face_gnum is simply referenced by d, so its lifecycle must be
     at least as long as that of d. */

  /* Face -> cell connectivity */
  /*---------------------------*/

  /* Build global cell numbering including parallel halos,
     except for periodic values */

  cell_gnum = cs_mesh_get_cell_gnum(mesh, 1);

  BFT_MALLOC(face_cell_g, n_faces*2, cs_gnum_t);

  for (i = 0; i < n_i_faces; i++) {
    face_cell_g[i*2] = cell_gnum[mesh->i_face_cells[i*2] - 1];
    face_cell_g[i*2 + 1] = cell_gnum[mesh->i_face_cells[i*2 + 1] - 1];
  }
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
    face_cell_g[j*2] = cell_gnum[mesh->b_face_cells[i] - 1];
    face_cell_g[j*2 + 1] = 0;
  }

  BFT_FREE(cell_gnum);

  /* Distribute to blocks and write */

  BFT_MALLOC(_face_cell_g,
             (face_bi.gnum_range[1] - face_bi.gnum_range[0]) * 2,
             cs_gnum_t);

  fvm_part_to_block_copy_array(d,
                               gnum_type,
                               2,
                               face_cell_g,
                               _face_cell_g);

  BFT_FREE(face_cell_g);

  cs_io_write_block_buffer("face_cells",
                           n_g_faces,
                           face_bi.gnum_range[0],
                           face_bi.gnum_range[1],
                           2, /* location_id, */
                           0, /* index id */
                           2, /* n_location_vals */
                           gnum_type,
                           _face_cell_g,
                           pp_out);

  BFT_FREE(_face_cell_g);

  /* Distribute and write blocks for group classes */
  /*-----------------------------------------------*/

  /* Now write pre-distributed blocks for cell group classes */

  cs_io_write_block_buffer("cell_group_class_id",
                           mesh->n_g_cells,
                           cell_bi.gnum_range[0],
                           cell_bi.gnum_range[1],
                           1, /* location_id, */
                           0, /* index id */
                           1, /* n_location_vals */
                           lnum_type,
                           _cell_gc_id,
                           pp_out);

  BFT_FREE(_cell_gc_id);

  /* Face group classes */

  BFT_MALLOC(face_gc_id, n_faces, cs_lnum_t);
  BFT_MALLOC(_face_gc_id,
             (face_bi.gnum_range[1] - face_bi.gnum_range[0]),
             cs_lnum_t);

  for (i = 0; i < n_i_faces; i++)
    face_gc_id[i] = mesh->i_face_family[i];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    face_gc_id[j] = mesh->b_face_family[i];

  /* Distribute to blocks and write */

  fvm_part_to_block_copy_array(d,
                               lnum_type,
                               1,
                               face_gc_id,
                               _face_gc_id);

  BFT_FREE(face_gc_id);

  cs_io_write_block_buffer("face_group_class_id",
                           n_g_faces,
                           face_bi.gnum_range[0],
                           face_bi.gnum_range[1],
                           2, /* location_id, */
                           0, /* index id */
                           1, /* n_location_vals */
                           lnum_type,
                           _face_gc_id,
                           pp_out);

  BFT_FREE(_face_gc_id);

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  _write_face_vertices_g(mesh, face_bi, d, pp_out);

  /* Free face part to block distribution structures */

  fvm_part_to_block_destroy(&d);
  BFT_FREE(face_gnum);

  /* Vertex coordinates */
  /*--------------------*/

  vtx_bi = fvm_part_to_block_compute_sizes(cs_glob_rank_id,
                                           cs_glob_n_ranks,
                                           0,
                                           0,
                                           mesh->n_g_vertices);

  d = fvm_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                       vtx_bi,
                                       mesh->n_vertices,
                                       mesh->global_vtx_num);

  BFT_MALLOC(_vtx_coords,
             (vtx_bi.gnum_range[1] - vtx_bi.gnum_range[0]) * 3,
             cs_real_t);

  fvm_part_to_block_copy_array(d,
                               real_type,
                               3,
                               mesh->vtx_coord,
                               _vtx_coords);

  cs_io_write_block_buffer("vertex_coords",
                           mesh->n_g_vertices,
                           vtx_bi.gnum_range[0],
                           vtx_bi.gnum_range[1],
                           3, /* location_id, */
                           0, /* index id */
                           3, /* n_location_vals */
                           real_type,
                           _vtx_coords,
                           pp_out);

  BFT_FREE(_vtx_coords);

  fvm_part_to_block_destroy(&d);
}

/*----------------------------------------------------------------------------
 * Write mesh periodicity data in parallel.
 *
 * parameters:
 *   perio_num       <-- periodicity number
 *   n_perio_couples <-- number of periodic face couples for this periodicity
 *   perio_couples   <-> periodic face couples for this periodicity
 *   pp_out          <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_perio_data_g(int         perio_num,
                         cs_lnum_t   n_perio_couples,
                         cs_gnum_t   perio_couples[],
                         cs_io_t    *pp_out)
{
  char section_name[32];
  fvm_part_to_block_info_t bi;

  cs_gnum_t  n_g_couples = 0;
  cs_gnum_t  *_perio_couples = NULL;
  const cs_gnum_t  *couple_g_num = NULL;
  fvm_part_to_block_t *d = NULL;

  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Create global couple numbering */

  fvm_io_num_t *c_io_num = fvm_io_num_create_from_adj_s(NULL,
                                                        perio_couples,
                                                        n_perio_couples,
                                                        2);

  n_g_couples = fvm_io_num_get_global_count(c_io_num);
  couple_g_num = fvm_io_num_get_global_num(c_io_num);

  /* Create associated block info and distribution */

  bi = fvm_part_to_block_compute_sizes(cs_glob_rank_id,
                                       cs_glob_n_ranks,
                                       0,
                                       0,
                                       n_g_couples);

  d = fvm_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                       bi,
                                       n_perio_couples,
                                       couple_g_num);

  BFT_MALLOC(_perio_couples,
             (bi.gnum_range[1] - bi.gnum_range[0]) * 2,
             cs_gnum_t);

  fvm_part_to_block_copy_array(d,
                               gnum_type,
                               2,
                               perio_couples,
                               _perio_couples);

  fvm_part_to_block_destroy(&d);

  c_io_num = fvm_io_num_destroy(c_io_num);

  /* Write face couples */

  sprintf(section_name, "periodicity_faces_%02d", perio_num);

  cs_io_write_block_buffer(section_name,
                           n_g_couples,
                           bi.gnum_range[0],
                           bi.gnum_range[1],
                           0, /* location_id, */
                           0, /* index id */
                           2, /* n_location_vals */
                           gnum_type,
                           _perio_couples,
                           pp_out);

  BFT_FREE(_perio_couples);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write mesh data in single-processor mode.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   pp_out   <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_data_l(const cs_mesh_t  *mesh,
                   cs_io_t          *pp_out)
{
  cs_lnum_t i, j, k, l, n_face_vertices, cell_id_0, cell_id_1;
  cs_gnum_t g_vtx_connect_size;

  cs_lnum_t *cell_gc_id = NULL, *face_gc_id = NULL;
  cs_gnum_t *face_vtx_idx_g = NULL, *face_vtx_g = NULL;
  cs_gnum_t *face_cell_g = NULL;
  cs_lnum_t *order = NULL, *i_order = NULL, *b_order = NULL;
  cs_real_t *vtx_coords = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t real_type
    = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t lnum_type
    = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;

  /* Prepare cell group classes */

  BFT_MALLOC(cell_gc_id, mesh->n_cells, cs_lnum_t);

  order = cs_order_gnum(NULL, mesh->global_cell_num, mesh->n_cells);
  for (i = 0; i < mesh->n_cells; i++)
    cell_gc_id[i] = mesh->cell_family[order[i]];
  BFT_FREE(order);

  /* Face ordering */

  i_order = cs_order_gnum(NULL, mesh->global_i_face_num, mesh->n_i_faces);
  if (mesh->n_b_faces > 0)
    b_order = cs_order_gnum(NULL, mesh->global_b_face_num, mesh->n_b_faces);

  /* Face -> cell connectivity (excluding periodic values )*/

  BFT_MALLOC(face_cell_g, n_faces*2, cs_gnum_t);

  if (mesh->global_cell_num != NULL) {

    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      cell_id_0 = mesh->i_face_cells[k*2] - 1;
      cell_id_1 = mesh->i_face_cells[k*2 + 1] - 1;
      if (cell_id_0 < mesh->n_cells)
        face_cell_g[i*2] = mesh->global_cell_num[cell_id_0];
      else
        face_cell_g[i*2] = 0;
      if (cell_id_1 < mesh->n_cells)
        face_cell_g[i*2 + 1] = mesh->global_cell_num[cell_id_1];
      else
        face_cell_g[i*2 + 1] = 0;
    }
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
      k = b_order[i];
      face_cell_g[j*2] = mesh->global_cell_num[mesh->b_face_cells[k] - 1];
      face_cell_g[j*2 + 1] = 0;
    }

  }
  else { /* if (mesh->global_cell_num == NULL) */

    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      cell_id_0 = mesh->i_face_cells[k*2] - 1;
      cell_id_1 = mesh->i_face_cells[k*2 + 1] - 1;
      if (cell_id_0 < mesh->n_cells)
        face_cell_g[i*2] = cell_id_0 + 1;
      else
        face_cell_g[i*2] = 0;
      if (cell_id_1 < mesh->n_cells)
        face_cell_g[i*2 + 1] = cell_id_1 + 1;
      else
        face_cell_g[i*2 + 1] = 0;
    }
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
      k = b_order[i];
      face_cell_g[j*2] = mesh->b_face_cells[k];
      face_cell_g[j*2 + 1] = 0;
    }

  }

  /* Write */

  cs_io_write_block_buffer("face_cells",
                           n_g_faces,
                           1,
                           n_g_faces + 1,
                           2, /* location_id, */
                           0, /* index id */
                           2, /* n_location_vals */
                           gnum_type,
                           face_cell_g,
                           pp_out);

  BFT_FREE(face_cell_g);

  /* Now write cell group classes */

  cs_io_write_block_buffer("cell_group_class_id",
                           mesh->n_g_cells,
                           1,
                           mesh->n_g_cells + 1,
                           1, /* location_id, */
                           0, /* index id */
                           1, /* n_location_vals */
                           lnum_type,
                           cell_gc_id,
                           pp_out);

  BFT_FREE(cell_gc_id);

  /* Face group classes */

  BFT_MALLOC(face_gc_id, n_faces, cs_lnum_t);

  for (i = 0; i < n_i_faces; i++)
    face_gc_id[i] = mesh->i_face_family[i_order[i]];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    face_gc_id[j] = mesh->b_face_family[b_order[i]];

  /* Write */

  cs_io_write_block_buffer("face_group_class_id",
                           n_g_faces,
                           1,
                           n_g_faces + 1,
                           2, /* location_id, */
                           0, /* index id */
                           1, /* n_location_vals */
                           lnum_type,
                           face_gc_id,
                           pp_out);

  BFT_FREE(face_gc_id);

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  BFT_MALLOC(face_vtx_idx_g, n_faces + 1, cs_gnum_t);

  face_vtx_idx_g[0] = 1;
  for (i = 0; i < n_i_faces; i++) {
    k = i_order[i];
    n_face_vertices = mesh->i_face_vtx_idx[k+1] - mesh->i_face_vtx_idx[k];
    face_vtx_idx_g[i+1] = face_vtx_idx_g[i] + n_face_vertices;
  }
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
    k = b_order[i];
    n_face_vertices = mesh->b_face_vtx_idx[k+1] - mesh->b_face_vtx_idx[k];
    face_vtx_idx_g[j+1] = face_vtx_idx_g[j] + n_face_vertices;
  }

  /* Now write buffer */

  cs_io_write_block_buffer("face_vertices_index",
                           n_g_faces + 1,
                           1,
                           n_faces + 2,
                           2, /* location_id, */
                           1, /* index id */
                           1, /* n_location_vals */
                           gnum_type,
                           face_vtx_idx_g,
                           pp_out);

  BFT_FREE(face_vtx_idx_g);

  /* Build connectivity */

  g_vtx_connect_size = (  mesh->i_face_vtx_connect_size
                        + mesh->b_face_vtx_connect_size);

  BFT_MALLOC(face_vtx_g, g_vtx_connect_size, cs_gnum_t);

  if (mesh->global_vtx_num != NULL) {
    l = 0;
    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      for (j = mesh->i_face_vtx_idx[k]; j < mesh->i_face_vtx_idx[k+1]; j++)
        face_vtx_g[l++]
          = mesh->global_vtx_num[mesh->i_face_vtx_lst[j - 1] - 1];
    }
    for (i = 0; i < n_b_faces; i++) {
      k = b_order[i];
      for (j = mesh->b_face_vtx_idx[k]; j < mesh->b_face_vtx_idx[k+1]; j++)
        face_vtx_g[l++]
          = mesh->global_vtx_num[mesh->b_face_vtx_lst[j - 1] - 1];
    }
  }
  else { /* if (mesh->global_vtx_num == NULL) */
    l = 0;
    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      for (j = mesh->i_face_vtx_idx[k]; j < mesh->i_face_vtx_idx[k+1]; j++)
        face_vtx_g[l++] = mesh->i_face_vtx_lst[j - 1];
    }
    for (i = 0; i < n_b_faces; i++) {
      k = b_order[i];
      for (j = mesh->b_face_vtx_idx[k]; j < mesh->b_face_vtx_idx[k+1]; j++)
        face_vtx_g[l++] = mesh->b_face_vtx_lst[j - 1];
    }
  }

  BFT_FREE(b_order);
  BFT_FREE(i_order);

  cs_io_write_block_buffer("face_vertices",
                           g_vtx_connect_size,
                           1,
                           g_vtx_connect_size + 1,
                           2, /* location_id, */
                           1, /* index id */
                           1, /* n_location_vals */
                           gnum_type,
                           face_vtx_g,
                           pp_out);

  BFT_FREE(face_vtx_g);

  /* Vertex coordinates */
  /*--------------------*/

  BFT_MALLOC(vtx_coords, mesh->n_vertices*3, cs_real_t);

  order = cs_order_gnum(NULL, mesh->global_vtx_num, mesh->n_vertices);
  for (i = 0; i < mesh->n_vertices; i++) {
    j = order[i];
    for (k = 0; k <3; k++)
      vtx_coords[i*3 + k] = mesh->vtx_coord[j*3 + k];
  }
  BFT_FREE(order);

  cs_io_write_block_buffer("vertex_coords",
                           mesh->n_g_vertices,
                           1,
                           mesh->n_vertices + 1,
                           3, /* location_id, */
                           0, /* index id */
                           3, /* n_location_vals */
                           real_type,
                           vtx_coords,
                           pp_out);

  BFT_FREE(vtx_coords);
}

/*----------------------------------------------------------------------------
 * Write mesh periodicity data in single-processor mode.
 *
 * parameters:
 *   perio_num       <-- periodicity number
 *   n_perio_couples <-- number of periodic face couples for this periodicity
 *   perio_couples   <-> periodic face couples for this periodicity
 *   pp_out          <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_perio_data_l(int         perio_num,
                         cs_lnum_t   n_perio_couples,
                         cs_gnum_t   perio_couples[],
                         cs_io_t    *pp_out)
{
  char section_name[32];

  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Write face couples */

  sprintf(section_name, "periodicity_faces_%02d", perio_num);

  cs_io_write_block_buffer(section_name,
                           n_perio_couples,
                           1,
                           n_perio_couples + 1,
                           0, /* location_id, */
                           0, /* index id */
                           2, /* n_location_vals */
                           gnum_type,
                           perio_couples,
                           pp_out);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Save a mesh as preprocessor data.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   filename <-- file name
 *----------------------------------------------------------------------------*/

void
cs_mesh_save(const cs_mesh_t  *mesh,
             const char       *filename)
{
  cs_lnum_t i;

  long  echo = CS_IO_ECHO_OPEN_CLOSE;
  cs_gnum_t g_i_face_vertices_size = 0, g_b_face_vertices_size = 0;
  cs_gnum_t g_face_vertices_size = 0;

  cs_gnum_t  *n_g_perio_faces = NULL;
  cs_lnum_t  *n_perio_faces = NULL;
  cs_gnum_t  **perio_faces = NULL;
  cs_io_t  *pp_out = NULL;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t lnum_type
    = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;

  /* Precompute some sizes */

  cs_mesh_g_face_vertices_sizes(mesh,
                                &g_i_face_vertices_size,
                                &g_b_face_vertices_size);

  g_face_vertices_size = g_i_face_vertices_size + g_b_face_vertices_size;

  /* Get periodic faces information if required */

  cs_mesh_get_perio_faces(mesh, &n_perio_faces, &perio_faces);

    /* In parallel, each periodic couple returned by
       cs_mesh_get_perio_faces() should appear on one rank only,
       so the global number of couples is simply the sum over all
       ranks of the local counts. */

  if (mesh->n_init_perio > 0)
    BFT_MALLOC(n_g_perio_faces, mesh->n_init_perio, cs_gnum_t);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t *_n_l_perio_faces = NULL;

    BFT_MALLOC(_n_l_perio_faces, mesh->n_init_perio, cs_gnum_t);

    for (i = 0; i < mesh->n_init_perio; i++)
      _n_l_perio_faces[i] = n_perio_faces[i];

    MPI_Allreduce(_n_l_perio_faces, n_g_perio_faces, mesh->n_init_perio,
                  CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);

    BFT_FREE(_n_l_perio_faces);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    for (i = 0; i < mesh->n_init_perio; i++)
      n_g_perio_faces[i] = n_perio_faces[i];
  }

  /* Open file for output */

#if defined(HAVE_MPI)
  pp_out = cs_io_initialize(filename,
                            "Face-based mesh definition, R0",
                            CS_IO_MODE_WRITE,
                            cs_glob_io_hints,
                            echo,
                            cs_glob_mpi_comm);
#else
  pp_out = cs_io_initialize(filename,
                            "Face-based mesh definition, R0",
                            CS_IO_MODE_WRITE,
                            echo,
                            -1);
#endif

  /* Write headers */
  /*---------------*/

  cs_io_write_global("start_block:dimensions", 0, 0, 0, 0, CS_DATATYPE_NULL,
                     NULL, pp_out);

  cs_io_write_global("n_cells", 1, 1, 0, 1, gnum_type,
                     &(mesh->n_g_cells),  pp_out);

  cs_io_write_global("n_faces", 1, 2, 0, 1, gnum_type,
                     &n_g_faces, pp_out);

  cs_io_write_global("n_vertices", 1, 3, 0, 1, gnum_type,
                     &(mesh->n_g_vertices), pp_out);

  cs_io_write_global("face_vertices_size", 1, 0, 0, 1, gnum_type,
                     &g_face_vertices_size, pp_out);

  cs_io_write_global("n_group_classes", 1, 0, 0, 1, lnum_type,
                     &(mesh->n_families), pp_out);

  cs_io_write_global("n_group_class_props_max", 1, 0, 0, 1, lnum_type,
                     &(mesh->n_max_family_items), pp_out);

  if (mesh->n_groups > 0) {

    cs_io_write_global("n_groups", 1, 0, 0, 1, lnum_type,
                       &(mesh->n_groups), pp_out);

    cs_io_write_global("group_name_index",
                       mesh->n_groups + 1, 0, 0, 1, lnum_type,
                       mesh->group_idx, pp_out);

    cs_io_write_global("group_name",
                       mesh->group_idx[mesh->n_groups] - 1, 0, 0, 1, CS_CHAR,
                       mesh->group_lst, pp_out);
  }

  cs_io_write_global("group_class_properties",
                     (mesh->n_families * mesh->n_max_family_items), 0, 0, 1,
                     lnum_type,
                     mesh->family_item, pp_out);

  if (mesh->n_init_perio > 0) {

    cs_int_t  n_rot_perio = 0;

    /* Count rotation periodicities */

    for (i = 0; i < mesh->n_init_perio; i++) {
      if (   fvm_periodicity_get_type(mesh->periodicity, i*2)
          >= FVM_PERIODICITY_ROTATION)
        n_rot_perio += 1;
    }

    /* Output periodicity information */

    cs_io_write_global("n_periodic_directions", 1, 0, 0, 1, lnum_type,
                       &(mesh->n_init_perio), pp_out);

    cs_io_write_global("n_periodic_rotations", 1, 0, 0, 1, lnum_type,
                       &n_rot_perio, pp_out);

    for (i = 0; i < mesh->n_init_perio; i++)
      n_g_perio_faces[i] *= 2;

    cs_io_write_global("n_periodic_faces", mesh->n_init_perio, 0, 0, 1,
                       gnum_type,
                       n_g_perio_faces, pp_out);

    for (i = 0; i < mesh->n_init_perio; i++)
      n_g_perio_faces[i] /= 2;

  }

  cs_io_write_global("end_block:dimensions", 0, 0, 0, 0, CS_DATATYPE_NULL,
                     NULL, pp_out);

  /* Write data */
  /*------------*/

  cs_io_write_global("start_block:data", 0, 0, 0, 0, CS_DATATYPE_NULL,
                     NULL, pp_out);

  /* Main mesh data */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _write_mesh_data_g(mesh, pp_out);

#endif

  if (cs_glob_n_ranks == 1)
    _write_mesh_data_l(mesh, pp_out);

  /* Periodcity data */

  for (i = 0; i < mesh->n_init_perio; i++) {

    _write_mesh_perio_metadata(mesh, i+1, pp_out);

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      _write_mesh_perio_data_g(i+1,
                               n_perio_faces[i],
                               perio_faces[i],
                               pp_out);
#endif

    if (cs_glob_n_ranks == 1)
      _write_mesh_perio_data_l(i+1,
                               n_perio_faces[i],
                               perio_faces[i],
                               pp_out);

    BFT_FREE(perio_faces[i]);

  }

  /* Close file */

  cs_io_write_global("end_block:data", 0, 0, 0, 0, CS_DATATYPE_NULL,
                     NULL, pp_out);

  if (n_perio_faces != NULL) {
    BFT_FREE(perio_faces);
    BFT_FREE(n_perio_faces);
    BFT_FREE(n_g_perio_faces);
  }

  cs_io_finalize(&pp_out);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
