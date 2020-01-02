/*============================================================================
 * \file Define cs_mesh_builder_t fields from cs_mesh_t fields.
 *============================================================================*/

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

#include "cs_base.h"
#include "cs_interface.h"
#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_to_builder.h"

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Define mesh builder's face vertices arrays in parallel mode.
 *
 * parameters:
 *    mesh <-- pointer to mesh structure
 *    mb   <-> pointer to mesh builder
 *    d    <-> face part to block distribution helper
 *----------------------------------------------------------------------------*/

static void
_mesh_to_builder_face_vertices_g(const cs_mesh_t     *mesh,
                                 cs_mesh_builder_t   *mb,
                                 cs_part_to_block_t  *d)
{
  cs_lnum_t i, j, k;
  cs_gnum_t block_size, n_block_faces, n_face_vertices;

  cs_lnum_t *face_vtx_idx = NULL;
  cs_gnum_t *face_vtx_g = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  BFT_FREE(mb->face_vertices_idx);
  BFT_FREE(mb->face_vertices);

  BFT_MALLOC(mb->face_vertices_idx,
             (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]) + 1,
             cs_lnum_t);

  BFT_MALLOC(face_vtx_idx, n_faces + 1, cs_lnum_t);

  face_vtx_idx[0] = 0;
  for (i = 0; i < n_i_faces; i++) {
    n_face_vertices = mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i];
    face_vtx_idx[i+1] = face_vtx_idx[i] + n_face_vertices;
  }
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
    n_face_vertices = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
    face_vtx_idx[j+1] = face_vtx_idx[j] + n_face_vertices;
  }

  cs_part_to_block_copy_index(d, face_vtx_idx, mb->face_vertices_idx);

  /* Build connectivity */

  BFT_MALLOC(face_vtx_g,
             mesh->i_face_vtx_connect_size + mesh->b_face_vtx_connect_size,
             cs_gnum_t);

  k = 0;
  for (i = 0; i < n_i_faces; i++) {
    for (j = mesh->i_face_vtx_idx[i]; j < mesh->i_face_vtx_idx[i+1]; j++)
      face_vtx_g[k++]
        = mesh->global_vtx_num[mesh->i_face_vtx_lst[j]];
  }
  for (i = 0; i < n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1]; j++)
      face_vtx_g[k++]
        = mesh->global_vtx_num[mesh->b_face_vtx_lst[j]];
  }

  n_block_faces = (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]);
  block_size = mb->face_vertices_idx[n_block_faces];

  BFT_MALLOC(mb->face_vertices, block_size, cs_gnum_t);

  cs_part_to_block_copy_indexed(d,
                                gnum_type,
                                face_vtx_idx,
                                face_vtx_g,
                                mb->face_vertices_idx,
                                mb->face_vertices);

  BFT_FREE(face_vtx_g);
  BFT_FREE(face_vtx_idx);
}

/*----------------------------------------------------------------------------
 * Write face->vertices connectivity in parallel.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   mb       <-> mesh builder
 *   transfer <-- if true, mesh transferred to builder; if false, builder
 *                is a temporary copy
 *   pp_out   <-- pointer to output file
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_g(const cs_mesh_t    *mesh,
                       cs_mesh_builder_t  *mb,
                       bool                transfer,
                       cs_io_t            *pp_out)
{
  cs_lnum_t i;
  cs_gnum_t block_size, g_vtx_connect_size, n_block_faces, n_face_vertices;
  cs_gnum_t idx_range[4];

  cs_gnum_t *face_vtx_idx_g = NULL;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  /* Copy index from block to global values */

  n_block_faces = (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]);
  BFT_MALLOC(face_vtx_idx_g, n_block_faces+ 1, cs_gnum_t);

  block_size = mb->face_vertices_idx[n_block_faces];
  MPI_Scan(&block_size, face_vtx_idx_g, 1, CS_MPI_GNUM, MPI_SUM,
           cs_glob_mpi_comm);
  face_vtx_idx_g[0] -= block_size;
  face_vtx_idx_g[0] += 1;

  for (i = 0; i < (cs_lnum_t)n_block_faces; i++) {
    n_face_vertices = mb->face_vertices_idx[i+1] - mb->face_vertices_idx[i];
    face_vtx_idx_g[i+1] = face_vtx_idx_g[i] + n_face_vertices;
  }

  idx_range[0] = mb->face_bi.gnum_range[0];
  idx_range[1] = mb->face_bi.gnum_range[1];
  if (mb->face_bi.gnum_range[0] >= n_g_faces) {
    idx_range[0] += 1;
    idx_range[1] += 1;
  }
  else if (mb->face_bi.gnum_range[1] >= n_g_faces + 1)
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

  if (transfer == true)
    cs_io_write_block("face_vertices",
                      g_vtx_connect_size,
                      idx_range[2],
                      idx_range[3],
                      2, /* location_id, */
                      1, /* index id */
                      1, /* n_location_vals */
                      gnum_type,
                      mb->face_vertices,
                      pp_out);
  else
    cs_io_write_block_buffer("face_vertices",
                             g_vtx_connect_size,
                             idx_range[2],
                             idx_range[3],
                             2, /* location_id, */
                             1, /* index id */
                             1, /* n_location_vals */
                             gnum_type,
                             mb->face_vertices,
                             pp_out);
}

/*----------------------------------------------------------------------------
 * Transfer or save mesh data to builder in parallel.
 *
 * parameters:
 *   mesh      <-> pointer to mesh structure
 *   mb        <-> mesh builder
 *   transfer  <-- if true, data is transferred from mesh to builder;
 *                 if false, builder fields are only used as a temporary
 *                 arrays.
 *   pp_out    <-> optional output file, or NULL
 *----------------------------------------------------------------------------*/

static void
_mesh_to_builder_g(cs_mesh_t          *mesh,
                   cs_mesh_builder_t  *mb,
                   bool                transfer,
                   cs_io_t            *pp_out)
{
  cs_lnum_t i, j;

  cs_part_to_block_t *d = NULL;

  cs_lnum_t *face_gc_id = NULL;
  cs_gnum_t *cell_gnum = NULL, *face_gnum = NULL;
  cs_gnum_t *face_cell_g = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_datatype_t real_type
    = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t int_type
    = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  /* Distribute cell group class info to blocks */
  /*---------------------------------------------*/

  BFT_MALLOC(mb->cell_gc_id,
             (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]),
             int);

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      mb->cell_bi,
                                      mesh->n_cells,
                                      mesh->global_cell_num);

  cs_part_to_block_copy_array(d,
                              int_type,
                              1,
                              mesh->cell_family,
                              mb->cell_gc_id);

  cs_part_to_block_destroy(&d);

  if (transfer == true)
    BFT_FREE(mesh->cell_family);

  /* Build global face part to block distribution structures */
  /*---------------------------------------------------------*/

  BFT_MALLOC(face_gnum, n_faces, cs_gnum_t);

  for (i = 0; i < n_i_faces; i++)
    face_gnum[i] = mesh->global_i_face_num[i];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    face_gnum[j] = mesh->global_b_face_num[i] + mesh->n_g_i_faces;

  if (transfer == true) {
    BFT_FREE(mesh->global_i_face_num);
    BFT_FREE(mesh->global_b_face_num);
  }

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      mb->face_bi,
                                      n_faces,
                                      face_gnum);
  cs_part_to_block_transfer_gnum(d, face_gnum);
  face_gnum = NULL;

  /* Face -> cell connectivity */
  /*---------------------------*/

  /* Build global cell numbering including parallel halos,
     except for periodic values */

  cell_gnum = cs_mesh_get_cell_gnum(mesh, 1);

  BFT_MALLOC(face_cell_g, n_faces*2, cs_gnum_t);

  for (i = 0; i < n_i_faces; i++) {
    face_cell_g[i*2]     = cell_gnum[mesh->i_face_cells[i][0]];
    face_cell_g[i*2 + 1] = cell_gnum[mesh->i_face_cells[i][1]];
  }
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
    face_cell_g[j*2] = cell_gnum[mesh->b_face_cells[i]];
    face_cell_g[j*2 + 1] = 0;
  }

  BFT_FREE(cell_gnum);

  if (transfer == true) {
    BFT_FREE(mesh->global_cell_num);
    BFT_FREE(mesh->i_face_cells);
    BFT_FREE(mesh->b_face_cells);
  }

  /* Distribute to blocks and write */

  BFT_MALLOC(mb->face_cells,
             (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]) * 2,
             cs_gnum_t);

  cs_part_to_block_copy_array(d,
                              gnum_type,
                              2,
                              face_cell_g,
                              mb->face_cells);

  BFT_FREE(face_cell_g);

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("face_cells",
                        mb->n_g_faces,
                        mb->face_bi.gnum_range[0],
                        mb->face_bi.gnum_range[1],
                        2, /* location_id, */
                        0, /* index id */
                        2, /* n_location_vals */
                        gnum_type,
                        mb->face_cells,
                        pp_out);
    else
      cs_io_write_block_buffer("face_cells",
                               mb->n_g_faces,
                               mb->face_bi.gnum_range[0],
                               mb->face_bi.gnum_range[1],
                               2, /* location_id, */
                               0, /* index id */
                               2, /* n_location_vals */
                               gnum_type,
                               mb->face_cells,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->face_cells);

  /* Now write pre-distributed blocks for cell group classes if required */

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("cell_group_class_id",
                        mesh->n_g_cells,
                        mb->cell_bi.gnum_range[0],
                        mb->cell_bi.gnum_range[1],
                        1, /* location_id, */
                        0, /* index id */
                        1, /* n_location_vals */
                        int_type,
                        mb->cell_gc_id,
                        pp_out);
    else
      cs_io_write_block_buffer("cell_group_class_id",
                               mesh->n_g_cells,
                               mb->cell_bi.gnum_range[0],
                               mb->cell_bi.gnum_range[1],
                               1, /* location_id, */
                               0, /* index id */
                               1, /* n_location_vals */
                               int_type,
                               mb->cell_gc_id,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->cell_gc_id);

  /* Face group classes */
  /*--------------------*/

  BFT_MALLOC(mb->face_gc_id,
             (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]),
             int);

  BFT_MALLOC(face_gc_id, n_faces, int);

  for (i = 0; i < n_i_faces; i++)
    face_gc_id[i] = mesh->i_face_family[i];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    face_gc_id[j] = mesh->b_face_family[i];

  if (transfer == true) {
    BFT_FREE(mesh->i_face_family);
    BFT_FREE(mesh->b_face_family);
  }

  /* Distribute to blocks, write if required */

  cs_part_to_block_copy_array(d,
                              int_type,
                              1,
                              face_gc_id,
                              mb->face_gc_id);

  BFT_FREE(face_gc_id);

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("face_group_class_id",
                        mb->n_g_faces,
                        mb->face_bi.gnum_range[0],
                        mb->face_bi.gnum_range[1],
                        2, /* location_id, */
                        0, /* index id */
                        1, /* n_location_vals */
                        int_type,
                        mb->face_gc_id,
                        pp_out);
    else
      cs_io_write_block_buffer("face_group_class_id",
                               mb->n_g_faces,
                               mb->face_bi.gnum_range[0],
                               mb->face_bi.gnum_range[1],
                               2, /* location_id, */
                               0, /* index id */
                               1, /* n_location_vals */
                               int_type,
                               mb->face_gc_id,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->face_gc_id);

  /* Face refinement generation */
  /*----------------------------*/

  if (mb->have_face_r_gen) {

    char *face_r_gen = NULL;

    BFT_MALLOC(mb->face_r_gen,
               (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]),
               char);

    BFT_MALLOC(face_r_gen, n_faces, char);

    for (i = 0; i < n_i_faces; i++)
      face_r_gen[i] = mesh->i_face_r_gen[i];
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
      face_r_gen[j] = 0;

    /* Distribute to blocks, write if required */

    cs_part_to_block_copy_array(d,
                                CS_CHAR,
                                1,
                                face_r_gen,
                                mb->face_r_gen);

    BFT_FREE(face_r_gen);

    if (pp_out != NULL) {
      if (transfer == true)
        cs_io_write_block("face_refinement_generation",
                          mb->n_g_faces,
                          mb->face_bi.gnum_range[0],
                          mb->face_bi.gnum_range[1],
                          2, /* location_id, */
                          0, /* index id */
                          1, /* n_location_vals */
                          CS_CHAR,
                          mb->face_r_gen,
                          pp_out);
      else
        cs_io_write_block_buffer("face_refinement_generation",
                                 mb->n_g_faces,
                                 mb->face_bi.gnum_range[0],
                                 mb->face_bi.gnum_range[1],
                                 2, /* location_id, */
                                 0, /* index id */
                                 1, /* n_location_vals */
                                 CS_CHAR,
                                 mb->face_r_gen,
                                 pp_out);

    }
  }

  if (transfer == true)
    BFT_FREE(mesh->i_face_r_gen);

  if (transfer == false)
    BFT_FREE(mb->face_r_gen);

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  _mesh_to_builder_face_vertices_g(mesh, mb, d);

  if (transfer == true) {
    BFT_FREE(mesh->i_face_vtx_idx);
    BFT_FREE(mesh->i_face_vtx_lst);
    BFT_FREE(mesh->b_face_vtx_idx);
    BFT_FREE(mesh->b_face_vtx_lst);
  }

  if (pp_out != NULL)
    _write_face_vertices_g(mesh, mb, transfer, pp_out);

  if (transfer == false) {
    BFT_FREE(mb->face_vertices_idx);
    BFT_FREE(mb->face_vertices);
  }

  /* Free face part to block distribution structures */

  cs_part_to_block_destroy(&d);

  /* Vertex coordinates */
  /*--------------------*/

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      mb->vertex_bi,
                                      mesh->n_vertices,
                                      mesh->global_vtx_num);

  BFT_MALLOC(mb->vertex_coords,
             (mb->vertex_bi.gnum_range[1] - mb->vertex_bi.gnum_range[0]) * 3,
             cs_real_t);

  cs_part_to_block_copy_array(d,
                              real_type,
                              3,
                              mesh->vtx_coord,
                              mb->vertex_coords);

  if (transfer == true)
    BFT_FREE(mesh->vtx_coord);

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("vertex_coords",
                        mesh->n_g_vertices,
                        mb->vertex_bi.gnum_range[0],
                        mb->vertex_bi.gnum_range[1],
                        3, /* location_id, */
                        0, /* index id */
                        3, /* n_location_vals */
                        real_type,
                        mb->vertex_coords,
                        pp_out);
    else
      cs_io_write_block_buffer("vertex_coords",
                               mesh->n_g_vertices,
                               mb->vertex_bi.gnum_range[0],
                               mb->vertex_bi.gnum_range[1],
                               3, /* location_id, */
                               0, /* index id */
                               3, /* n_location_vals */
                               real_type,
                               mb->vertex_coords,
                               pp_out);
  }

  if (transfer == true)
    BFT_FREE(mesh->global_vtx_num);
  else
    BFT_FREE(mb->vertex_coords);

  if (transfer == true)
    BFT_FREE(mesh->global_vtx_num);
  else
    BFT_FREE(mb->vertex_coords);

  cs_part_to_block_destroy(&d);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Transfer or save mesh data to builder in serial mode.
 *
 * parameters:
 *   mesh      <-> pointer to mesh structure
 *   mb        <-> mesh builder
 *   transfer  <-- if true, data is transferred from mesh to builder;
 *                 if false, builder fields are only used as a temporary
 *                 arrays.
 *   pp_out    <-> optional output file, or NULL
 *----------------------------------------------------------------------------*/

static void
_mesh_to_builder_l(cs_mesh_t          *mesh,
                   cs_mesh_builder_t  *mb,
                   bool                transfer,
                   cs_io_t            *pp_out)
{
  cs_lnum_t i, j, k, l, n_face_vertices, cell_id_0, cell_id_1;
  cs_gnum_t g_vtx_connect_size;

  cs_lnum_t *order = NULL, *i_order = NULL, *b_order = NULL;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_faces = mesh->n_i_faces + mesh->n_b_faces;

  const cs_datatype_t real_type
    = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t int_type
    = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  /* Face ordering */

  i_order = cs_order_gnum(NULL, mesh->global_i_face_num, mesh->n_i_faces);
  if (mesh->n_b_faces > 0)
    b_order = cs_order_gnum(NULL, mesh->global_b_face_num, mesh->n_b_faces);

  if (transfer == true) {
    BFT_FREE(mesh->global_i_face_num);
    BFT_FREE(mesh->global_b_face_num);
  }

  /* Face -> cell connectivity (excluding periodic values )*/

  BFT_MALLOC(mb->face_cells, n_faces*2, cs_gnum_t);

  if (mesh->global_cell_num != NULL) {

    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      cell_id_0 = mesh->i_face_cells[k][0];
      cell_id_1 = mesh->i_face_cells[k][1];
      if (cell_id_0 < mesh->n_cells)
        mb->face_cells[i*2] = mesh->global_cell_num[cell_id_0];
      else
        mb->face_cells[i*2] = 0;
      if (cell_id_1 < mesh->n_cells)
        mb->face_cells[i*2 + 1] = mesh->global_cell_num[cell_id_1];
      else
        mb->face_cells[i*2 + 1] = 0;
    }
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
      k = b_order[i];
      mb->face_cells[j*2] = mesh->global_cell_num[mesh->b_face_cells[k]];
      mb->face_cells[j*2 + 1] = 0;
    }

    if (transfer == true)
      BFT_FREE(mesh->global_cell_num);

  }
  else { /* if (mesh->global_cell_num == NULL) */

    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      cell_id_0 = mesh->i_face_cells[k][0];
      cell_id_1 = mesh->i_face_cells[k][1];
      if (cell_id_0 < mesh->n_cells)
        mb->face_cells[i*2] = cell_id_0 + 1;
      else
        mb->face_cells[i*2] = 0;
      if (cell_id_1 < mesh->n_cells)
        mb->face_cells[i*2 + 1] = cell_id_1 + 1;
      else
        mb->face_cells[i*2 + 1] = 0;
    }
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
      k = b_order[i];
      mb->face_cells[j*2] = mesh->b_face_cells[k] + 1;
      mb->face_cells[j*2 + 1] = 0;
    }

  }

  if (transfer == true) {
    BFT_FREE(mesh->i_face_cells);
    BFT_FREE(mesh->b_face_cells);
  }

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("face_cells",
                        mb->n_g_faces,
                        1,
                        mb->n_g_faces + 1,
                        2, /* location_id, */
                        0, /* index id */
                        2, /* n_location_vals */
                        gnum_type,
                        mb->face_cells,
                        pp_out);
    else
      cs_io_write_block_buffer("face_cells",
                               mb->n_g_faces,
                               1,
                               mb->n_g_faces + 1,
                               2, /* location_id, */
                               0, /* index id */
                               2, /* n_location_vals */
                               gnum_type,
                               mb->face_cells,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->face_cells);

  /* Cell group classes */

  BFT_MALLOC(mb->cell_gc_id, mesh->n_cells, int);

  order = cs_order_gnum(NULL, mesh->global_cell_num, mesh->n_cells);
  for (i = 0; i < mesh->n_cells; i++)
    mb->cell_gc_id[i] = mesh->cell_family[order[i]];
  BFT_FREE(order);

  if (transfer == true)
    BFT_FREE(mesh->cell_family);

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("cell_group_class_id",
                        mesh->n_g_cells,
                        1,
                        mesh->n_g_cells + 1,
                        1, /* location_id, */
                        0, /* index id */
                        1, /* n_location_vals */
                        int_type,
                        mb->cell_gc_id,
                        pp_out);
    else
      cs_io_write_block_buffer("cell_group_class_id",
                               mesh->n_g_cells,
                               1,
                               mesh->n_g_cells + 1,
                               1, /* location_id, */
                               0, /* index id */
                               1, /* n_location_vals */
                               int_type,
                               mb->cell_gc_id,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->cell_gc_id);

  /* Face group classes */

  BFT_MALLOC(mb->face_gc_id, n_faces, int);

  for (i = 0; i < n_i_faces; i++)
    mb->face_gc_id[i] = mesh->i_face_family[i_order[i]];
  for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
    mb->face_gc_id[j] = mesh->b_face_family[b_order[i]];

  if (transfer == true) {
    BFT_FREE(mesh->i_face_family);
    BFT_FREE(mesh->b_face_family);
  }

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("face_group_class_id",
                        mb->n_g_faces,
                        1,
                        mb->n_g_faces + 1,
                        2, /* location_id, */
                        0, /* index id */
                        1, /* n_location_vals */
                        int_type,
                        mb->face_gc_id,
                        pp_out);
    else
      cs_io_write_block_buffer("face_group_class_id",
                               mb->n_g_faces,
                               1,
                               mb->n_g_faces + 1,
                               2, /* location_id, */
                               0, /* index id */
                               1, /* n_location_vals */
                               int_type,
                               mb->face_gc_id,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->face_gc_id);

  /* Face refinement generation */

  if (mb->have_face_r_gen) {

    BFT_MALLOC(mb->face_r_gen, n_faces, char);

    for (i = 0; i < n_i_faces; i++)
      mb->face_r_gen[i] = mesh->i_face_family[i_order[i]];
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++)
      mb->face_r_gen[j] = mesh->b_face_family[b_order[i]];

    if (transfer == true)
      BFT_FREE(mesh->i_face_r_gen);

    if (pp_out != NULL) {
      if (transfer == true)
        cs_io_write_block("face_refinement_generation",
                          mb->n_g_faces,
                          1,
                          mb->n_g_faces + 1,
                          2, /* location_id, */
                          0, /* index id */
                          1, /* n_location_vals */
                          CS_CHAR,
                          mb->face_r_gen,
                          pp_out);
      else
        cs_io_write_block_buffer("face_refinement_generation",
                                 mb->n_g_faces,
                                 1,
                                 mb->n_g_faces + 1,
                                 2, /* location_id, */
                                 0, /* index id */
                                 1, /* n_location_vals */
                                 CS_CHAR,
                                 mb->face_r_gen,
                                 pp_out);
    }

    if (transfer == false)
      BFT_FREE(mb->face_r_gen);
  }

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  if (pp_out != NULL) {

    cs_gnum_t *face_vtx_idx_g = NULL;

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
                             mb->n_g_faces + 1,
                             1,
                             n_faces + 2,
                             2, /* location_id, */
                             1, /* index id */
                             1, /* n_location_vals */
                             gnum_type,
                             face_vtx_idx_g,
                             pp_out);

    BFT_FREE(face_vtx_idx_g);

  }

  if (transfer == true) {

    BFT_MALLOC(mb->face_vertices_idx, n_faces + 1, cs_lnum_t);

    mb->face_vertices_idx[0] = 0;
    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      n_face_vertices = mesh->i_face_vtx_idx[k+1] - mesh->i_face_vtx_idx[k];
      mb->face_vertices_idx[i+1] = mb->face_vertices_idx[i] + n_face_vertices;
    }
    for (i = 0, j = n_i_faces; i < n_b_faces; i++, j++) {
      k = b_order[i];
      n_face_vertices = mesh->b_face_vtx_idx[k+1] - mesh->b_face_vtx_idx[k];
      mb->face_vertices_idx[j+1] = mb->face_vertices_idx[j] + n_face_vertices;
    }

  }

  /* Build connectivity */

  g_vtx_connect_size = (  mesh->i_face_vtx_connect_size
                        + mesh->b_face_vtx_connect_size);

  BFT_MALLOC(mb->face_vertices, g_vtx_connect_size, cs_gnum_t);

  if (mesh->global_vtx_num != NULL) {
    l = 0;
    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      for (j = mesh->i_face_vtx_idx[k]; j < mesh->i_face_vtx_idx[k+1]; j++)
        mb->face_vertices[l++]
          = mesh->global_vtx_num[mesh->i_face_vtx_lst[j - 1] - 1];
    }
    for (i = 0; i < n_b_faces; i++) {
      k = b_order[i];
      for (j = mesh->b_face_vtx_idx[k]; j < mesh->b_face_vtx_idx[k+1]; j++)
        mb->face_vertices[l++]
          = mesh->global_vtx_num[mesh->b_face_vtx_lst[j]];
    }
  }
  else { /* if (mesh->global_vtx_num == NULL) */
    l = 0;
    for (i = 0; i < n_i_faces; i++) {
      k = i_order[i];
      for (j = mesh->i_face_vtx_idx[k]; j < mesh->i_face_vtx_idx[k+1]; j++)
        mb->face_vertices[l++] = mesh->i_face_vtx_lst[j] + 1;
    }
    for (i = 0; i < n_b_faces; i++) {
      k = b_order[i];
      for (j = mesh->b_face_vtx_idx[k]; j < mesh->b_face_vtx_idx[k+1]; j++)
        mb->face_vertices[l++] = mesh->b_face_vtx_lst[j] + 1;
    }
  }

  BFT_FREE(b_order);
  BFT_FREE(i_order);

  if (transfer == true) {
    BFT_FREE(mesh->i_face_vtx_idx);
    BFT_FREE(mesh->i_face_vtx_lst);
    BFT_FREE(mesh->b_face_vtx_idx);
    BFT_FREE(mesh->b_face_vtx_lst);
  }

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("face_vertices",
                        g_vtx_connect_size,
                        1,
                        g_vtx_connect_size + 1,
                        2, /* location_id, */
                        1, /* index id */
                        1, /* n_location_vals */
                        gnum_type,
                        mb->face_vertices,
                        pp_out);
    else
      cs_io_write_block_buffer("face_vertices",
                               g_vtx_connect_size,
                               1,
                               g_vtx_connect_size + 1,
                               2, /* location_id, */
                               1, /* index id */
                               1, /* n_location_vals */
                               gnum_type,
                               mb->face_vertices,
                               pp_out);
  }

  if (transfer == false)
    BFT_FREE(mb->face_vertices);

  /* Vertex coordinates */
  /*--------------------*/

  BFT_MALLOC(mb->vertex_coords, mesh->n_vertices*3, cs_real_t);

  order = cs_order_gnum(NULL, mesh->global_vtx_num, mesh->n_vertices);
  for (i = 0; i < mesh->n_vertices; i++) {
    j = order[i];
    for (k = 0; k <3; k++)
      mb->vertex_coords[i*3 + k] = mesh->vtx_coord[j*3 + k];
  }
  BFT_FREE(order);

  if (transfer == true)
    BFT_FREE(mesh->vtx_coord);

  if (pp_out != NULL) {
    if (transfer == true)
      cs_io_write_block("vertex_coords",
                        mesh->n_g_vertices,
                        1,
                        mesh->n_vertices + 1,
                        3, /* location_id, */
                        0, /* index id */
                        3, /* n_location_vals */
                        real_type,
                        mb->vertex_coords,
                        pp_out);
    else
      cs_io_write_block_buffer("vertex_coords",
                               mesh->n_g_vertices,
                               1,
                               mesh->n_vertices + 1,
                               3, /* location_id, */
                               0, /* index id */
                               3, /* n_location_vals */
                               real_type,
                               mb->vertex_coords,
                               pp_out);
  }

  if (transfer == true)
    BFT_FREE(mesh->global_vtx_num);
  else
    BFT_FREE(mb->vertex_coords);
}

/*----------------------------------------------------------------------------
 * Save a mesh as preprocessor data.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   mb       <-- pointer to optional mesh builder structure, or NULL
 *   pp_out   <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_dimensions(cs_mesh_t          *mesh,
                  cs_mesh_builder_t  *mb,
                  cs_io_t            *pp_out)
{
  cs_lnum_t i;

  const cs_gnum_t n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;
  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  const cs_datatype_t lnum_type
    = (sizeof(cs_lnum_t) == 8) ? CS_INT64 : CS_INT32;

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
                     &(mb->n_g_face_connect_size), pp_out);

  cs_io_write_global("n_group_classes", 1, 0, 0, 1, lnum_type,
                     &(mesh->n_families), pp_out);

  cs_io_write_global("n_group_class_props_max", 1, 0, 0, 1, lnum_type,
                     &(mesh->n_max_family_items), pp_out);

  if (mesh->n_groups > 0) {

    cs_io_write_global("n_groups", 1, 0, 0, 1, lnum_type,
                       &(mesh->n_groups), pp_out);

    cs_lnum_t *_group_idx;
    BFT_MALLOC(_group_idx, mesh->n_groups + 1, cs_lnum_t);
    for (i = 0; i < mesh->n_groups + 1; i++)
      _group_idx[i] = mesh->group_idx[i] + 1;
    cs_io_write_global("group_name_index",
                       mesh->n_groups + 1, 0, 0, 1, lnum_type,
                       _group_idx, pp_out);
    BFT_FREE(_group_idx);

    cs_io_write_global("group_name",
                       mesh->group_idx[mesh->n_groups], 0, 0, 1, CS_CHAR,
                       mesh->group, pp_out);
  }

  cs_io_write_global("group_class_properties",
                     (mesh->n_families * mesh->n_max_family_items), 0, 0, 1,
                     lnum_type,
                     mesh->family_item, pp_out);

  if (mesh->n_init_perio > 0) {

    int n_rot_perio = 0;

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
      mb->n_g_per_face_couples[i] *= 2;

    cs_io_write_global("n_periodic_faces", mesh->n_init_perio, 0, 0, 1,
                       gnum_type,
                       mb->n_g_per_face_couples, pp_out);

    for (i = 0; i < mesh->n_init_perio; i++)
      mb->n_g_per_face_couples[i] /= 2;

  }

  cs_io_write_global("end_block:dimensions", 0, 0, 0, 0, CS_DATATYPE_NULL,
                     NULL, pp_out);
}

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
 * Write mesh periodicity data in parallel.
 *
 * parameters:
 *   perio_num       <-- periodicity number
 *   n_perio_couples <-- number of periodic face couples for this periodicity
 *   perio_couples   <-> periodic face couples for this periodicity
 *   min_rank_step   <-- minimum rank step between blocks
 *   transfer        <-- if true, mesh transferred to builder;
 *                       if false, builder is a temporary copy
 *   pp_out          <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_perio_data_g(int         perio_num,
                         cs_lnum_t   n_perio_couples,
                         cs_gnum_t   perio_couples[],
                         int         min_rank_step,
                         bool        transfer,
                         cs_io_t    *pp_out)
{
  char section_name[32];
  cs_block_dist_info_t bi;

  cs_gnum_t  n_g_couples = 0;
  cs_gnum_t  *_perio_couples = NULL;
  const cs_gnum_t  *couple_g_num = NULL;
  cs_part_to_block_t *d = NULL;

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

  bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                   cs_glob_n_ranks,
                                   min_rank_step,
                                   0,
                                   n_g_couples);

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      bi,
                                      n_perio_couples,
                                      couple_g_num);

  BFT_MALLOC(_perio_couples,
             (bi.gnum_range[1] - bi.gnum_range[0]) * 2,
             cs_gnum_t);

  cs_part_to_block_copy_array(d,
                              gnum_type,
                              2,
                              perio_couples,
                              _perio_couples);

  cs_part_to_block_destroy(&d);

  c_io_num = fvm_io_num_destroy(c_io_num);

  /* Write face couples */

  sprintf(section_name, "periodicity_faces_%02d", perio_num);

  if (transfer == true)
    cs_io_write_block(section_name,
                      n_g_couples,
                      bi.gnum_range[0],
                      bi.gnum_range[1],
                      0, /* location_id, */
                      0, /* index id */
                      2, /* n_location_vals */
                      gnum_type,
                      _perio_couples,
                      pp_out);
  else
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
 * Write mesh periodicity data in single-processor mode.
 *
 * parameters:
 *   perio_num       <-- periodicity number
 *   n_perio_couples <-- number of periodic face couples for this periodicity
 *   perio_couples   <-> periodic face couples for this periodicity
 *   transfer        <-- if true, mesh transferred to builder;
 *                       if false, builder is a temporary copy
 *   pp_out          <-> output file
 *----------------------------------------------------------------------------*/

static void
_write_mesh_perio_data_l(int         perio_num,
                         cs_lnum_t   n_perio_couples,
                         cs_gnum_t   perio_couples[],
                         bool        transfer,
                         cs_io_t    *pp_out)
{
  char section_name[32];

  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Write face couples */

  sprintf(section_name, "periodicity_faces_%02d", perio_num);

  if (transfer == true)
    cs_io_write_block(section_name,
                      n_perio_couples,
                      1,
                      n_perio_couples + 1,
                      0, /* location_id, */
                      0, /* index id */
                      2, /* n_location_vals */
                      gnum_type,
                      perio_couples,
                      pp_out);
  else
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer mesh to mesh builder structure.
 *
 * As the dataflow is very similar, but may be done array-by array to minimize
 * memory overhead, this function also handles a part of the output
 * to file needed to save a mesh file.
 *
 * \param[in, out]  mesh      pointer to mesh structure
 * \param[in, out]  mb        pointer to mesh builder structure
 * \param[in]       transfer  if true, data is transferred from mesh to builder;
 *                            if false, builder fields are only used as a
 *                            temporary arrays.
 * \param[in, out]  pp_out    optional output file, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_to_builder(cs_mesh_t          *mesh,
                   cs_mesh_builder_t  *mb,
                   bool                transfer,
                   cs_io_t            *pp_out)
{
  int i;
  cs_gnum_t g_i_face_vertices_size = 0, g_b_face_vertices_size = 0;

  /* Clear rebuildable mesh data in case of transfer;
     halos are kept at this stage, as they are required for this algorithm */

  if (transfer)
    cs_mesh_free_rebuildable(mesh, false);

  /* Clear previous builder data if present (periodicity done separately) */

  BFT_FREE(mb->face_cells);
  BFT_FREE(mb->face_vertices_idx);
  BFT_FREE(mb->face_vertices);
  BFT_FREE(mb->cell_gc_id);
  BFT_FREE(mb->face_gc_id);
  BFT_FREE(mb->face_r_gen);
  BFT_FREE(mb->vertex_coords);

  /* Precompute some sizes */

  mb->n_g_faces = mesh->n_g_i_faces + mesh->n_g_b_faces;

  cs_mesh_builder_define_block_dist(mb,
                                    cs_glob_rank_id,
                                    cs_glob_n_ranks,
                                    mb->min_rank_step,
                                    0,
                                    mesh->n_g_cells,
                                    mb->n_g_faces,
                                    mesh->n_g_vertices);

  cs_mesh_g_face_vertices_sizes(mesh,
                                &g_i_face_vertices_size,
                                &g_b_face_vertices_size);

  mb->n_g_face_connect_size = g_i_face_vertices_size + g_b_face_vertices_size;

  /* Get refinement info if needed */

  int r_flag = 0;
  for (cs_lnum_t j = 0; j < mesh->n_i_faces; j++) {
    if (mesh->i_face_r_gen[j] != 0) {
      r_flag = 1;
      break;
    }
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    int _r_flag = r_flag;
    MPI_Allreduce(&_r_flag, &r_flag, 1,
                  MPI_INT, MPI_MAX, cs_glob_mpi_comm);
  }
#endif

  mb->have_face_r_gen = (r_flag) ? true : false;

  /* Get periodic faces information if required */

  cs_mesh_to_builder_perio_faces(mesh, mb);

  /* Write metadata if output is required */

  if (pp_out != NULL)
    _write_dimensions(mesh, mb, pp_out);

  /* Main mesh data */

  if (pp_out != NULL)
    cs_io_write_global("start_block:data", 0, 0, 0, 0, CS_DATATYPE_NULL,
                       NULL, pp_out);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _mesh_to_builder_g(mesh, mb, transfer, pp_out);

#endif

  if (cs_glob_n_ranks == 1)
    _mesh_to_builder_l(mesh, mb, transfer, pp_out);

  /* Finish clearing rebuildable mesh data in case of transfer  */

  if (transfer) {
    if (mesh->vtx_interfaces != NULL)
      cs_interface_set_destroy(&(mesh->vtx_interfaces));
    if (mesh->halo != NULL)
      cs_halo_destroy(&(mesh->halo));
  }

  /* Periodicity data */

  if (pp_out != NULL) {

    for (i = 0; i < mesh->n_init_perio; i++) {

      _write_mesh_perio_metadata(mesh, i+1, pp_out);

#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1)
        _write_mesh_perio_data_g(i+1,
                                 mb->n_per_face_couples[i],
                                 mb->per_face_couples[i],
                                 mb->min_rank_step,
                                 transfer,
                                 pp_out);
#endif

      if (cs_glob_n_ranks == 1)
        _write_mesh_perio_data_l(i+1,
                                 mb->n_per_face_couples[i],
                                 mb->per_face_couples[i],
                                 transfer,
                                 pp_out);

      if (transfer == false)
        BFT_FREE(mb->per_face_couples[i]);

    }

    cs_io_write_global("end_block:data", 0, 0, 0, 0, CS_DATATYPE_NULL,
                       NULL, pp_out);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer mesh partitioning info to mesh builder structure.
 *
 * \param[in]       mesh      pointer to mesh structure
 * \param[in, out]  mb        pointer to mesh builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_to_builder_partition(const cs_mesh_t    *mesh,
                             cs_mesh_builder_t  *mb)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    const cs_datatype_t int_type
      = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

    /* Distribute cell group class info to blocks */
    /*---------------------------------------------*/

    mb->cell_bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                              cs_glob_n_ranks,
                                              mb->min_rank_step,
                                              0,
                                              mesh->n_g_cells);

    mb->have_cell_rank = true;
    BFT_REALLOC(mb->cell_rank,
                (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]),
                int);

    int *cell_rank;
    BFT_MALLOC(cell_rank, mesh->n_cells, int);
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      cell_rank[i] = cs_glob_rank_id;

    cs_part_to_block_t *d
      = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                        mb->cell_bi,
                                        mesh->n_cells,
                                        mesh->global_cell_num);

    cs_part_to_block_copy_array(d,
                                int_type,
                                1,
                                cell_rank,
                                mb->cell_rank);

    cs_part_to_block_destroy(&d);

    BFT_FREE(cell_rank);

  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct periodic faces info from mesh to builder.
 *
 * \param[in]       mesh   pointer to mesh structure
 * \param[in, out]  mb     pointer to mesh builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_to_builder_perio_faces(const cs_mesh_t    *mesh,
                               cs_mesh_builder_t  *mb)
{
  cs_lnum_t i;

  /* Get periodic faces information if required */

  mb->n_perio = mesh->n_init_perio;

  if (mesh->n_init_perio < 1)
    return;

  cs_mesh_get_perio_faces(mesh,
                          &(mb->n_per_face_couples),
                          &(mb->per_face_couples));

  /* In parallel, each periodic couple returned by
     cs_mesh_get_perio_faces() should appear on one rank only,
     so the global number of couples is simply the sum over all
     ranks of the local counts. */

  BFT_MALLOC(mb->n_g_per_face_couples, mesh->n_init_perio, cs_gnum_t);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t *_n_l_perio_faces = NULL;

    BFT_MALLOC(_n_l_perio_faces, mesh->n_init_perio, cs_gnum_t);

    for (i = 0; i < mesh->n_init_perio; i++)
      _n_l_perio_faces[i] = mb->n_per_face_couples[i];

    MPI_Allreduce(_n_l_perio_faces, mb->n_g_per_face_couples, mesh->n_init_perio,
                  CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);

    BFT_FREE(_n_l_perio_faces);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    for (i = 0; i < mesh->n_init_perio; i++)
      mb->n_g_per_face_couples[i] = mb->n_per_face_couples[i];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
