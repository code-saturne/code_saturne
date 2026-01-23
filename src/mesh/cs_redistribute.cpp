/*============================================================================
 * Repartitioning and redistribution mesh data and fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <ctime>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_all_to_all.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_field.h"
#include "base/cs_mem.h"
#include "base/cs_io.h"
#include "base/cs_order.h"
#include "base/cs_renumber.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_volume_zone.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_builder.h"
#include "mesh/cs_mesh_coherency.h"
#include "mesh/cs_mesh_from_builder.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_to_builder.h"

#include "alge/cs_cell_to_vertex.h"
#include "alge/cs_gradient.h"
#include "alge/cs_matrix_default.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_redistribute.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare strings (bsearch function).
 *
 * parameters:
 *   a <-> pointer to first number
 *   b <-> pointer to second number
 *
 * returns:
 *   -1 if a < b, 0 if a == b, 1 if a > b
 *----------------------------------------------------------------------------*/

static int
_cmp_gnum(const void  *a,
          const void  *b)
{
  const cs_gnum_t *_a = (const cs_gnum_t *)a;
  const cs_gnum_t *_b = (const cs_gnum_t *)b;
  if (*_a < *_b) return -1;
  if (*_a > *_b) return 1;
  return 0;
}

/*----------------------------------------------------------------------------
 * Determine if a given face should be distributed.
 *
 * parameters:
 *   mesh      <-- pointer to mesh
 *   f_id      <-- face id
 *   rank      <-- local rank
 *   cell_rank <-- adjacent cell rank
 *
 * returns:
 *   true if face should be redistributed.
 *----------------------------------------------------------------------------*/

static bool
_distribute_face(cs_mesh_t *mesh,
                 cs_lnum_t f_id,
                 int rank,
                 int cell_rank[])
{
  bool distribute = true;

  for (int j = 0; j < 2; j++) {
    cs_lnum_t c_id = mesh->i_face_cells[f_id][j];

    // proc face ?
    if (c_id >= mesh->n_cells) {

      // do not distribute a proc face if remote proc id is lower than local one
      if (cell_rank[c_id] < rank) {
        distribute = false;
      }

      break;
    }
  }

  return distribute;
}

/*----------------------------------------------------------------------------
 * Exchange and order integer arrays.
 *
 * parameters:
 *   d            <-- pointer to parallel distributor
 *   stride       <-- face id
 *   send_data    <-- values to send
 *   recv_data    <-- values to receive
 *   b_face_order <-- boundary faces ordering by associated global number
 *   b_face_n2o   <-- boundary faces new to old mapping
 *
 * returns:
 *   true if face should be redistributed.
 *----------------------------------------------------------------------------*/

template <typename T>
static void
_exchange_and_order_b(cs_all_to_all_t  *d,
                      int               stride,
                      bool              reverse,
                      const T           send_data[],
                      T                 recv_data[],
                      const cs_lnum_t   b_face_order[],
                      const cs_lnum_t   b_face_n2o[])
{
  cs_lnum_t n_elem_recv = cs_all_to_all_n_elts_dest(d);

  cs_all_to_all_copy_array(d,
                           stride,
                           reverse,
                           send_data,
                           recv_data);

  if (!b_face_order && !b_face_n2o) return;

  T *copy;
  CS_MALLOC(copy, n_elem_recv*stride, T);
  memcpy(copy, recv_data, n_elem_recv*stride*sizeof(T));

  for (cs_lnum_t new_slot = 0; new_slot < n_elem_recv; new_slot++) {
    cs_lnum_t old_slot = b_face_n2o ? b_face_n2o[new_slot] : new_slot;
    cs_lnum_t pre_slot = b_face_order ? b_face_order[old_slot] : old_slot;
    const T *orig = copy + stride*pre_slot;
    T *dest = recv_data + stride*new_slot;
    for (int i = 0; i < stride; i++)
      dest[i] = orig[i];
  }

  CS_FREE(copy);
}

/*----------------------------------------------------------------------------
 * Distribute BC coefficients.
 *
 * parameters:
 *   bfd          <-- pointer to parallel distributor
 *   n_elem       <-- number of elements
 *   stride       <-- face id
 *   copy_buf     <-- copy buffer
 *   coeff        <-> BC coefficients
 *   b_face_order <-- boundary faces ordering by associated global number
 *   b_face_n2o   <-- boundary faces new to old mapping
 *
 * returns:
 *   true if face should be redistributed.
 *----------------------------------------------------------------------------*/

static void
_distribute_bc_coeff(cs_all_to_all_t  *bfd,
                     cs_lnum_t         n_elem,
                     int               stride,
                     cs_real_t         copy_buf[],
                     cs_real_t        *coeff[],
                     const cs_lnum_t   b_face_order[],
                     const cs_lnum_t   b_face_n2o[])
{
  if (*coeff) {
    memcpy(copy_buf, *coeff, n_elem*stride*sizeof(cs_real_t));
    CS_FREE(*coeff);
  }

  cs_lnum_t n_elts = cs_all_to_all_n_elts_dest(bfd);

  CS_REALLOC(*coeff, n_elts*stride, cs_real_t);

  _exchange_and_order_b(bfd,
                        stride,
                        false,
                        copy_buf,
                        *coeff,
                        b_face_order,
                        b_face_n2o);
}

/*----------------------------------------------------------------------------
 * Distribute BC types.
 *
 * parameters:
 *   bfd           <-- pointer to parallel distributor
 *   n_b_faces_ini <-- initial number of boundary faces
 *   b_face_order  <-- boundary faces ordering by associated global number
 *   b_face_n2o    <-- boundary faces new to old mapping
 *
 * returns:
 *   true if face should be redistributed.
 *----------------------------------------------------------------------------*/

static void
_distribute_bc_type(cs_all_to_all_t  *bfd,
                    cs_lnum_t         n_b_faces_ini,
                    const cs_lnum_t   b_face_order[],
                    const cs_lnum_t   b_face_n2o[])
{
  cs_mesh_t *mesh = cs_glob_mesh;
  cs_boundary_condition_pm_info_t *pm = cs_glob_bc_pm_info;

  int *copy;
  CS_MALLOC(copy, n_b_faces_ini, int);

  memcpy(copy, pm->izfppp, n_b_faces_ini*sizeof(int));
  CS_REALLOC(pm->izfppp, mesh->n_b_faces, int);
  _exchange_and_order_b(bfd, 1, false, copy, pm->izfppp,
                        b_face_order, b_face_n2o);

  if (pm->iautom) {
    memcpy(copy, pm->iautom, n_b_faces_ini*sizeof(int));
    CS_REALLOC(pm->iautom, mesh->n_b_faces, int);
    _exchange_and_order_b(bfd, 1, false, copy, pm->iautom,
                          b_face_order, b_face_n2o);
  }

  int **bc_type = nullptr;
  cs_boundary_conditions_get_bc_type_addr(&bc_type);
  memcpy(copy, *bc_type, n_b_faces_ini*sizeof(int));
  CS_REALLOC(*bc_type, mesh->n_b_faces, int);
  _exchange_and_order_b(bfd, 1, false, copy, *bc_type,
                        b_face_order, b_face_n2o);

  CS_FREE(copy);

  cs_glob_bc_type = *bc_type;
}

/*----------------------------------------------------------------------------
 * Determine destination ranks for cells based on refinment info
 *
 * parameters:
 *   mesh           <-- pointer to mesh
 *   cell_dest_rank --> cell destination rank
 *----------------------------------------------------------------------------*/

void
_compute_cell_dest_rank(cs_mesh_t  *mesh,
                        int         cell_dest[])
{
  // Cells sharing a parent should have the same destination.

  char *c_r_level;
  CS_MALLOC(c_r_level, mesh->n_cells, char);
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
    c_r_level[c_id] = 0;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    // Can't merge an original face.
    if (mesh->i_face_r_gen[f_id] < 1)
      continue;

    for (int i = 0; i < 2; i++) {
      cs_lnum_t c_id = mesh->i_face_cells[f_id][i];
      if (c_id < mesh->n_cells) {
        // cell refinement level is max i_face_r_gen over its internal faces.
        if (mesh->i_face_r_gen[f_id] > c_r_level[c_id]) {
          c_r_level[c_id] = mesh->i_face_r_gen[f_id];
        }
      }
    }
  }

  int reloop = 0;

  cs_lnum_t *c_equiv = nullptr;
  CS_MALLOC(c_equiv, mesh->n_cells, cs_lnum_t);
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
    c_equiv[c_id] = c_id;

  do {

    reloop = 0;

    for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
      if (mesh->i_face_r_gen[f_id] < 1)
        continue;

      cs_lnum_t c_id0 = mesh->i_face_cells[f_id][0];
      cs_lnum_t c_id1 = mesh->i_face_cells[f_id][1];

      if (c_id0 >= mesh->n_cells || c_id1 >= mesh->n_cells)
        continue;

      if (mesh->i_face_r_gen[f_id] == c_r_level[c_id0] &&
          mesh->i_face_r_gen[f_id] == c_r_level[c_id1]) {
        cs_lnum_t min_equiv = cs::min(c_equiv[c_id0], c_equiv[c_id1]);

        if (c_equiv[c_id0] > min_equiv) {
          c_equiv[c_id0] = min_equiv;
          reloop = 1;
        }

        if (c_equiv[c_id1] > min_equiv) {
          c_equiv[c_id1] = min_equiv;
          reloop = 1;
        }
      }
    }
  } while (reloop);

  CS_FREE(c_r_level);

  for (int c_id = 0; c_id < mesh->n_cells; c_id++) {
    if (c_id == c_equiv[c_id]) {
      cell_dest[c_id] = rand() % cs_glob_n_ranks;
    } else {
      assert(c_id > c_equiv[c_id]);
      cell_dest[c_id] = cell_dest[c_equiv[c_id]];
    }
  }

  CS_FREE(c_equiv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Redistribute mesh and fields based on a cell destination rank map.
 *
 * If no cell map is given, a random one is created.
 *
 * \param[in]  cell_dest_rank  destination rank for each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_redistribute(int  *cell_dest_rank)
{
  static int iter = 0;
  iter++;

  cs_mesh_t *mesh = cs_glob_mesh;

  MPI_Comm comm = cs_glob_mpi_comm;

  const cs_lnum_t n_cells_ini = mesh->n_cells;
  const cs_lnum_t n_b_faces_ini = mesh->n_b_faces;

  // create block-to-part builder
  bool transfer = false;
  cs_io_t *output = nullptr;
  cs_mesh_builder_t *builder = cs_mesh_builder_create();

  // create the mesh builder.
  // Note: the builder is still useful throughout the part-to-part distribution
  // for the following steps:
  // - cs_partition uses it to compute the block-to-part cell map,
  // - block-to-part vertex distribution (for now),
  // - periodicity info for halo construction once the mesh is redistributed.

  cs_mesh_to_builder(mesh, builder, transfer, output);

  /* Cells
     ----- */

  cs_all_to_all_t *cd = nullptr;
  int cd_flags = 0;

  // build the cell part-to-part distributor

  cs_lnum_t n_cells = 0;

  // Random cell destinations if cell_dest_rank is not prescribed.

  int *_dest_rank = nullptr;
  if (cell_dest_rank == nullptr) {

    CS_MALLOC(_dest_rank, mesh->n_cells_with_ghosts, int);
    cell_dest_rank = _dest_rank;

    // RNG
    static bool rng_on = false;
    if (!rng_on) {
      rng_on = true;

      unsigned int seed;

      if (cs_glob_rank_id <= 0) {
        seed = time(NULL);
        //seed = 1; // Fix the seed for debugging.
      }

      MPI_Bcast(&seed, 1, MPI_UINT32_T, 0, comm);

      srand(seed);
    }

    cs_lnum_t g_min_n_cells;

    do {
      if (cd) cs_all_to_all_destroy(&cd);

      _compute_cell_dest_rank(mesh, cell_dest_rank);

      cd = cs_all_to_all_create(mesh->n_cells,
                                cd_flags,
                                nullptr,
                                cell_dest_rank,
                                comm);
      n_cells = cs_all_to_all_n_elts_dest(cd);

      MPI_Allreduce(&n_cells, &g_min_n_cells, 1, MPI_INT, MPI_MIN, comm);

    } while (g_min_n_cells == 0);
  }
  else {
    cd = cs_all_to_all_create(mesh->n_cells,
                              cd_flags,
                              nullptr,
                              cell_dest_rank,
                              comm);
    n_cells = cs_all_to_all_n_elts_dest(cd);
  }

  // distribute global cell numbers

  cs_gnum_t *cell_gnum = cs_all_to_all_copy_array(cd,
                                                  1,
                                                  false,
                                                  mesh->global_cell_num);

  // Sort global cell numbers in increasing order
  // this order will be used for the exchange of any cell-centered data.

  cs_lnum_t *cell_order = cs_order_gnum(nullptr, cell_gnum, n_cells);

  cs_gnum_t *global_cell_num;
  CS_MALLOC(global_cell_num, n_cells, cs_gnum_t);
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    global_cell_num[i] = cell_gnum[cell_order[i]];
  }

  CS_FREE(cell_gnum);

  // Distribute cell family

  int *cell_family = cs_all_to_all_copy_array(cd, 1, false, mesh->cell_family);
  CS_REALLOC(mesh->cell_family, n_cells, int);
  for (cs_lnum_t i = 0; i < n_cells; i++)
    mesh->cell_family[i] = cell_family[cell_order[i]];
  CS_FREE(cell_family);

  /* Boundary faces
     -------------- */

  // Boundary face destination/
  int *b_face_dest;
  CS_MALLOC(b_face_dest, mesh->n_b_faces, int);
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    b_face_dest[f_id] = cell_dest_rank[mesh->b_face_cells[f_id]];

  // Boundary face distributor.
  int bfd_flags = 0;
  cs_all_to_all_t *bfd = cs_all_to_all_create(mesh->n_b_faces,
                                              bfd_flags,
                                              nullptr,
                                              b_face_dest,
                                              comm);

  cs_all_to_all_transfer_dest_rank(bfd, &b_face_dest);

  // New number of boundary faces.
  cs_lnum_t n_b_faces = cs_all_to_all_n_elts_dest(bfd);

  // Global boundary face numbers distribution.
  cs_gnum_t *b_face_gnum;
  b_face_gnum = cs_all_to_all_copy_array(bfd,
                                         1,
                                         false,
                                         mesh->global_b_face_num);

  cs_lnum_t *b_face_order = cs_order_gnum(nullptr, b_face_gnum, n_b_faces);

  CS_REALLOC(mesh->global_b_face_num, n_b_faces, cs_gnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    mesh->global_b_face_num[f_id] = b_face_gnum[b_face_order[f_id]];
  CS_FREE(b_face_gnum);

  // AMR fields

  if (mesh->b_face_r_c_idx) {
    int *b_face_r_c_idx_s = nullptr;
    CS_MALLOC(b_face_r_c_idx_s, mesh->n_b_faces, int);
    for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
      b_face_r_c_idx_s[f_id] = mesh->b_face_r_c_idx[f_id];
    CS_FREE(mesh->b_face_r_c_idx);
    int *b_face_r_c_idx_r = cs_all_to_all_copy_array(bfd,
                                                     1,
                                                     false,
                                                     b_face_r_c_idx_s);
    CS_FREE(b_face_r_c_idx_s);
    if (b_face_r_c_idx_r) {
      CS_MALLOC(mesh->b_face_r_c_idx, n_b_faces, char);
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
        mesh->b_face_r_c_idx[f_id] = b_face_r_c_idx_r[b_face_order[f_id]];
      CS_FREE(b_face_r_c_idx_r);
    }
  }

  // Boundary face family.

  int *b_face_family = cs_all_to_all_copy_array(bfd,
                                                1,
                                                false,
                                                mesh->b_face_family);
  CS_REALLOC(mesh->b_face_family, n_b_faces, int);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    mesh->b_face_family[f_id] = b_face_family[b_face_order[f_id]];
  CS_FREE(b_face_family);

  // Boundary face cells distribution.

  cs_gnum_t *g_b_face_cells_s;
  CS_MALLOC(g_b_face_cells_s, mesh->n_b_faces, cs_gnum_t);
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    g_b_face_cells_s[f_id] = mesh->global_cell_num[mesh->b_face_cells[f_id]];
  CS_FREE(mesh->b_face_cells);

  cs_gnum_t *g_b_face_cells_r = cs_all_to_all_copy_array(bfd,
                                                         1,
                                                         false,
                                                         g_b_face_cells_s);
  CS_FREE(g_b_face_cells_s);

  // Transform boundary face-cell connectivity into local indices.

  cs_lnum_t *b_face_cells;
  CS_MALLOC(b_face_cells, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_gnum_t gc = g_b_face_cells_r[b_face_order[f_id]];
    void *found = bsearch(&gc, global_cell_num, n_cells, sizeof(cs_gnum_t),
      _cmp_gnum);
    assert(found);
    b_face_cells[f_id] = (cs_gnum_t *)found - global_cell_num;
  }
  CS_FREE(g_b_face_cells_r);
  mesh->b_face_cells = b_face_cells;

  // Boundary face-vertex connectivity
  // for now, vertices are global.

  for (cs_lnum_t f_id = mesh->n_b_faces; f_id >= 1; f_id--)
    mesh->b_face_vtx_idx[f_id] -= mesh->b_face_vtx_idx[f_id-1];

  cs_lnum_t *b_face_vtx_idx_r;
  CS_MALLOC(b_face_vtx_idx_r, n_b_faces+1, cs_lnum_t);

  cs_all_to_all_copy_array(bfd,
                           1,
                           false,
                           mesh->b_face_vtx_idx+1,
                           b_face_vtx_idx_r+1);

  // Apply b_face_order on the indices.

  cs_lnum_t *b_face_vtx_idx;
  CS_MALLOC(b_face_vtx_idx, n_b_faces+1, cs_lnum_t);
  b_face_vtx_idx[0] = 0;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    b_face_vtx_idx[f_id+1] = b_face_vtx_idx_r[b_face_order[f_id]+1];
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    b_face_vtx_idx[f_id+1] += b_face_vtx_idx[f_id];

  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    mesh->b_face_vtx_idx[f_id+1] += mesh->b_face_vtx_idx[f_id];

  b_face_vtx_idx_r[0] = 0;
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    b_face_vtx_idx_r[f_id+1] += b_face_vtx_idx_r[f_id];

  assert(b_face_vtx_idx_r[n_b_faces] == b_face_vtx_idx[n_b_faces]);

  cs_gnum_t *g_b_face_vtx_lst_s;
  CS_MALLOC(g_b_face_vtx_lst_s, mesh->b_face_vtx_connect_size, cs_gnum_t);
  for (cs_lnum_t i = 0; i < mesh->b_face_vtx_connect_size; i++)
    g_b_face_vtx_lst_s[i] = mesh->global_vtx_num[mesh->b_face_vtx_lst[i]];
  CS_FREE(mesh->b_face_vtx_lst);

  cs_gnum_t *g_b_face_vtx_lst_r = cs_all_to_all_copy_indexed(bfd,
                                                             false,
                                                             mesh->b_face_vtx_idx,
                                                             g_b_face_vtx_lst_s,
                                                             b_face_vtx_idx_r);

  CS_FREE(mesh->b_face_vtx_idx);
  CS_FREE(g_b_face_vtx_lst_s);

  /* Interior faces
     -------------- */

  // We need the destination of the halo cells.
  cs_halo_sync_untyped(mesh->halo, mesh->halo_type, sizeof(int), cell_dest_rank);

  // We need the ranks of the halo cells.
  int *cell_rank;
  CS_MALLOC(cell_rank, mesh->n_cells_with_ghosts, int);
  for (int i = 0; i < mesh->n_cells; i++) cell_rank[i] = cs_glob_rank_id;
  cs_halo_sync_untyped(mesh->halo, mesh->halo_type, sizeof(int), cell_rank);

  // List of internal faces to distribute, and their destination.
  // one internal face can appear twice in the list.
  cs_lnum_t *i_face_list;
  cs_gnum_t *i_face_gnum_s;
  int *i_face_dest;

  CS_MALLOC(i_face_list, 2*mesh->n_i_faces, cs_lnum_t);
  CS_MALLOC(i_face_gnum_s, 2*mesh->n_i_faces, cs_gnum_t);
  CS_MALLOC(i_face_dest, 2*mesh->n_i_faces, int);

  // Number of internal faces to send.
  cs_lnum_t n_i_f_s = 0;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {

    if (!_distribute_face(mesh, f_id, cs_glob_rank_id, cell_rank)) continue;

    int dest_rank = -1;
    int prev_rank = -1;

    for (int j = 0; j < 2; j++) {
      dest_rank = cell_dest_rank[mesh->i_face_cells[f_id][j]];
      if (dest_rank != prev_rank) {
        i_face_list[n_i_f_s] = f_id;
        i_face_dest[n_i_f_s] = dest_rank;
        i_face_gnum_s[n_i_f_s] = mesh->global_i_face_num[f_id];
        n_i_f_s++;
        prev_rank = dest_rank;
      }
    }
  }

  CS_FREE(mesh->global_i_face_num);
  CS_FREE(cell_rank);

  // Interior faces distributor.
  int ifd_flags = 0;
  cs_all_to_all_t *ifd = cs_all_to_all_create(n_i_f_s,
                                              ifd_flags,
                                              nullptr,
                                              i_face_dest,
                                              comm);

  cs_all_to_all_transfer_dest_rank(ifd, &i_face_dest);

  // New number of internal faces.
  cs_lnum_t n_i_faces = cs_all_to_all_n_elts_dest(ifd);

  // Exchange global internal face numbers.
  cs_gnum_t *i_face_gnum = cs_all_to_all_copy_array(ifd,
                                                    1,
                                                    false,
                                                    i_face_gnum_s);
  CS_FREE(i_face_gnum_s);

  cs_lnum_t *i_face_order = cs_order_gnum(nullptr, i_face_gnum, n_i_faces);

  CS_MALLOC(mesh->global_i_face_num, n_i_faces, cs_gnum_t);
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    mesh->global_i_face_num[i] = i_face_gnum[i_face_order[i]];
  }
  CS_FREE(i_face_gnum);

  // Internal face family.
  int *i_face_family_s;
  CS_MALLOC(i_face_family_s, n_i_f_s, int);
  for (cs_lnum_t i = 0; i < n_i_f_s; i++) {
    i_face_family_s[i] = mesh->i_face_family[i_face_list[i]];
  }
  CS_FREE(mesh->i_face_family);
  int *i_face_family_r = cs_all_to_all_copy_array(ifd,
                                                  1,
                                                  false,
                                                  i_face_family_s);
  CS_FREE(i_face_family_s);
  CS_MALLOC(mesh->i_face_family, n_i_faces, int);
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    mesh->i_face_family[i] = i_face_family_r[i_face_order[i]];
  }
  CS_FREE(i_face_family_r);

  // AMR fields
  int *i_face_r_gen_s = nullptr;

  if (mesh->i_face_r_gen) {
    CS_MALLOC(i_face_r_gen_s, n_i_f_s, int);
    for (cs_lnum_t i = 0; i < n_i_f_s; i++) {
      cs_lnum_t f_id = i_face_list[i];
      i_face_r_gen_s[i] = mesh->i_face_r_gen[f_id];
    }
    CS_FREE(mesh->i_face_r_gen);
  }

  int *i_face_r_gen_r = cs_all_to_all_copy_array(ifd,
                                                 1,
                                                 false,
                                                 i_face_r_gen_s);

  CS_FREE(i_face_r_gen_s);
  if (i_face_r_gen_r) {

    CS_MALLOC(mesh->i_face_r_gen, n_i_faces, char);
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
      mesh->i_face_r_gen[f_id] = i_face_r_gen_r[i_face_order[f_id]];

    CS_FREE(i_face_r_gen_r);
  }

  // Distribute the (global) internal face-cell connectivity.
  int blank_perio = 1;
  cell_gnum = cs_mesh_get_cell_gnum(mesh, blank_perio);

  cs_gnum_t *g_i_face_cells_s;
  CS_MALLOC(g_i_face_cells_s, 2*n_i_f_s, cs_gnum_t);

  for (cs_lnum_t i = 0; i < n_i_f_s; i++) {
    cs_lnum_t f_id = i_face_list[i];
    for (int j = 0; j < 2; j++)
      g_i_face_cells_s[2*i+j] = cell_gnum[mesh->i_face_cells[f_id][j]];
  }

  CS_FREE(cell_gnum);

  cs_gnum_t *g_i_face_cells_r = cs_all_to_all_copy_array(ifd,
                                                         2,
                                                         false,
                                                         g_i_face_cells_s);

  CS_FREE(g_i_face_cells_s);

  cs_gnum_t *g_i_face_cells;
  CS_MALLOC(g_i_face_cells, 2*n_i_faces, cs_gnum_t);
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    for (int j = 0; j < 2; j++)
      g_i_face_cells[2*i+j] = g_i_face_cells_r[2*i_face_order[i]+j];
  }
  CS_FREE(g_i_face_cells_r);

  // Convert face-cells connectivity into local indices.

  cs_lnum_2_t *i_face_cells;
  CS_MALLOC(i_face_cells, n_i_faces, cs_lnum_2_t);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    for (int i = 0; i < 2; i++) {
      cs_gnum_t gc = g_i_face_cells[2*f_id+i];
      void *found = bsearch(&gc, global_cell_num, n_cells, sizeof(cs_gnum_t),
        _cmp_gnum);
      i_face_cells[f_id][i] = found ?
                              (cs_gnum_t *)found - global_cell_num :
                              -1;
    }
  }
  CS_FREE(g_i_face_cells);
  CS_FREE(mesh->i_face_cells);
  mesh->i_face_cells = i_face_cells;

  // Face-vertex connectivity index.
  cs_lnum_t *i_face_vtx_idx_s;
  CS_MALLOC(i_face_vtx_idx_s, n_i_f_s+1, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_i_f_s; i++) {
    cs_lnum_t f_id = i_face_list[i];
    i_face_vtx_idx_s[i+1] = mesh->i_face_vtx_idx[f_id+1] -
                            mesh->i_face_vtx_idx[f_id];
  }

  cs_lnum_t *i_face_vtx_idx_r;
  CS_MALLOC(i_face_vtx_idx_r, n_i_faces+1, cs_lnum_t);
  cs_all_to_all_copy_array(ifd,
                           1,
                           false,
                           i_face_vtx_idx_s+1,
                           i_face_vtx_idx_r+1);

  // Apply i_face_order on the indices.
  cs_lnum_t *i_face_vtx_idx;
  CS_MALLOC(i_face_vtx_idx, n_i_faces+1, cs_lnum_t);
  i_face_vtx_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_i_faces; i++)
    i_face_vtx_idx[i+1] = i_face_vtx_idx_r[i_face_order[i]+1];
  for (cs_lnum_t i = 0; i < n_i_faces; i++)
    i_face_vtx_idx[i+1] += i_face_vtx_idx[i];

  i_face_vtx_idx_s[0] = 0;
  for (cs_lnum_t i = 0; i < n_i_f_s; i++)
    i_face_vtx_idx_s[i+1] += i_face_vtx_idx_s[i];

  i_face_vtx_idx_r[0] = 0;
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    i_face_vtx_idx_r[f_id+1] += i_face_vtx_idx_r[f_id];

  assert(i_face_vtx_idx[n_i_faces] == i_face_vtx_idx_r[n_i_faces]);

  // Face-vertex connectivity list;
  // for now, exchange the global vertex indices.
  cs_gnum_t *g_i_face_vtx_lst_s;
  CS_MALLOC(g_i_face_vtx_lst_s, i_face_vtx_idx_s[n_i_f_s], cs_gnum_t);
  for (cs_lnum_t i = 0; i < n_i_f_s; i++) {
    cs_gnum_t *ptr = g_i_face_vtx_lst_s + i_face_vtx_idx_s[i];
    cs_lnum_t f_id = i_face_list[i];
    for (cs_lnum_t j = mesh->i_face_vtx_idx[f_id];
         j < mesh->i_face_vtx_idx[f_id+1];
         j++) {
      *ptr++ = mesh->global_vtx_num[mesh->i_face_vtx_lst[j]];
    }
  }
  CS_FREE(mesh->i_face_vtx_lst);
  CS_FREE(mesh->global_vtx_num);

  cs_gnum_t *g_i_face_vtx_lst_r = cs_all_to_all_copy_indexed(ifd,
                                                             false,
                                                             i_face_vtx_idx_s,
                                                             g_i_face_vtx_lst_s,
                                                             i_face_vtx_idx_r);

  CS_FREE(i_face_vtx_idx_s);
  CS_FREE(g_i_face_vtx_lst_s);
  CS_FREE(mesh->i_face_vtx_idx);

  /* Vertices
     -------- */

  size_t n_vertices = 0;
  cs_gnum_t *global_vtx_num = nullptr;

  cs_gnum_t *all_vtx;
  CS_MALLOC(all_vtx, i_face_vtx_idx_r[n_i_faces]+b_face_vtx_idx_r[n_b_faces],
            cs_gnum_t);
  for (cs_lnum_t i = 0; i < i_face_vtx_idx_r[n_i_faces]; i++)
    all_vtx[i] = g_i_face_vtx_lst_r[i];
  for (cs_lnum_t i = 0; i < b_face_vtx_idx_r[n_b_faces]; i++)
    all_vtx[i+i_face_vtx_idx_r[n_i_faces]] = g_b_face_vtx_lst_r[i];

  cs_order_single_gnum(i_face_vtx_idx_r[n_i_faces]+b_face_vtx_idx_r[n_b_faces],
                       1,
                       all_vtx,
                       &n_vertices,
                       &global_vtx_num);
  CS_FREE(all_vtx);

  mesh->global_vtx_num = global_vtx_num;

  // Note: to simplify implementation for now, we use the block vertices.
  // Later on, we might do a "pure" part-to-part distribution for vertices,
  // for example by splitting the vertices into two part: inner vertices
  // (shared by pure internal/boundary faces only), and outer vertices
  // (those that are connected to faces shared by neighbouring procs).
  // TODO(Imad): distribute vertex-centered fields.

  cs_all_to_all_t *vd
    = cs_all_to_all_create_from_block(n_vertices,
                                      CS_ALL_TO_ALL_USE_DEST_ID,
                                      global_vtx_num,
                                      builder->vertex_bi,
                                      comm);

  CS_FREE(mesh->vtx_coord);
  mesh->vtx_coord = cs_all_to_all_copy_array(vd,
                                             3,
                                             true,
                                             builder->vertex_coords);

  CS_FREE(builder->vertex_coords);

  if (mesh->have_r_gen) {
    assert(builder->vtx_r_gen);
    mesh->vtx_r_gen = cs_all_to_all_copy_array(vd,
                                               1,
                                               true,
                                               builder->vtx_r_gen);
  }
  CS_FREE(builder->vtx_r_gen);

  // Copy vertex-centered fields here.

  cs_all_to_all_destroy(&vd);

  mesh->n_vertices = n_vertices;

  // Apply face order on the received face-vtx connectivities.

  cs_gnum_t *g_b_face_vtx_lst;
  CS_MALLOC(g_b_face_vtx_lst, b_face_vtx_idx[n_b_faces], cs_gnum_t);
  memset(g_b_face_vtx_lst, 0, b_face_vtx_idx[n_b_faces]*sizeof(cs_gnum_t));

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_gnum_t *ptr = g_b_face_vtx_lst + b_face_vtx_idx[i];
    for (cs_lnum_t j = b_face_vtx_idx_r[b_face_order[i]];
         j < b_face_vtx_idx_r[b_face_order[i]+1];
         j++) {
      *ptr++ = g_b_face_vtx_lst_r[j];
    }
  }
  CS_FREE(g_b_face_vtx_lst_r);
  CS_FREE(b_face_vtx_idx_r);

  cs_gnum_t *g_i_face_vtx_lst;
  CS_MALLOC(g_i_face_vtx_lst, i_face_vtx_idx[n_i_faces], cs_gnum_t);
  memset(g_i_face_vtx_lst, 0, i_face_vtx_idx[n_i_faces]*sizeof(cs_gnum_t));

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    cs_gnum_t *ptr = g_i_face_vtx_lst + i_face_vtx_idx[i];
    for (cs_lnum_t j = i_face_vtx_idx_r[i_face_order[i]];
         j < i_face_vtx_idx_r[i_face_order[i]+1]; j++) {
      *ptr++ = g_i_face_vtx_lst_r[j];
    }
  }
  CS_FREE(g_i_face_vtx_lst_r);
  CS_FREE(i_face_vtx_idx_r);

  // Convert face-vertex connectivity into local indices.

  cs_lnum_t *i_face_vtx_lst;
  CS_MALLOC(i_face_vtx_lst, i_face_vtx_idx[n_i_faces], cs_lnum_t);
  for (cs_lnum_t i = 0; i < i_face_vtx_idx[n_i_faces]; i++) {
    cs_gnum_t gv = g_i_face_vtx_lst[i];
    // TODO: test inlined code_saturne function such as _g_id_binary_find
    // instead of bsearch to reduce overhead of function calls.
    void *found = bsearch(&gv, global_vtx_num, n_vertices, sizeof(cs_gnum_t),
                          _cmp_gnum);
    assert(found);
    i_face_vtx_lst[i] = (cs_gnum_t *)found - global_vtx_num;
  }
  CS_FREE(g_i_face_vtx_lst);
  mesh->n_i_faces = n_i_faces;
  mesh->i_face_vtx_idx = i_face_vtx_idx;
  mesh->i_face_vtx_lst = i_face_vtx_lst;
  mesh->i_face_vtx_connect_size = mesh->i_face_vtx_idx[mesh->n_i_faces];

  cs_lnum_t *b_face_vtx_lst;
  CS_MALLOC(b_face_vtx_lst, b_face_vtx_idx[n_b_faces], cs_lnum_t);
  for (cs_lnum_t i = 0; i < b_face_vtx_idx[n_b_faces]; i++) {
    cs_gnum_t gv = g_b_face_vtx_lst[i];
    void *found = bsearch(&gv, global_vtx_num, n_vertices, sizeof(cs_gnum_t),
      _cmp_gnum);
    assert(found);
    b_face_vtx_lst[i] = (cs_gnum_t *)found - global_vtx_num;
  }
  CS_FREE(g_b_face_vtx_lst);
  mesh->n_b_faces = n_b_faces;
  mesh->b_face_vtx_idx = b_face_vtx_idx;
  mesh->b_face_vtx_lst = b_face_vtx_lst;
  mesh->b_face_vtx_connect_size = mesh->b_face_vtx_idx[mesh->n_b_faces];
  mesh->n_b_faces_all = mesh->n_b_faces;

  mesh->n_cells = n_cells;
  CS_FREE(mesh->global_cell_num);
  mesh->global_cell_num = global_cell_num;

  // Note: attributes not updated by re-partition (conf. cs_mesh.h):
  // - n_g_i_c_faces
  // - periodicity features
  // - face status flags

  // update the halo
  cs_mesh_builder_destroy(&builder);
  cs_mesh_free_rebuildable(mesh, true);
  mesh->n_cells_with_ghosts = mesh->n_cells;
  cs_mesh_init_halo(mesh, nullptr, mesh->halo_type, -1, true);
  cs_mesh_update_auxiliary(mesh);

  cs_lnum_t *cell_n2o, *i_face_n2o, *b_face_n2o, *vtx_n2o;
  cell_n2o = i_face_n2o = b_face_n2o = vtx_n2o = nullptr;
  // cs_renumber_mesh(mesh, &cell_n2o, &i_face_n2o, &b_face_n2o, &vtx_n2o);

  {
    if (!cell_n2o) {
      mesh->cell_numbering = cs_numbering_create_default(mesh->n_cells);
      mesh->i_face_numbering = cs_numbering_create_default(mesh->n_i_faces);
      mesh->b_face_numbering = cs_numbering_create_default(mesh->n_b_faces);
      mesh->vtx_numbering = cs_numbering_create_default(mesh->n_vertices);
    }
  }

  cs_mesh_init_group_classes(mesh);

  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);
  cs_mesh_quantities_compute(mesh, cs_glob_mesh_quantities);

  cs_mesh_update_selectors(mesh);
  cs_mesh_location_build(mesh, -1);
  bool mesh_modified = true;
  cs_volume_zone_build_all(mesh_modified);
  cs_boundary_zone_build_all(mesh_modified);

  cs_ext_neighborhood_reduce(mesh, cs_glob_mesh_quantities);

  _distribute_bc_type(bfd,
                      n_b_faces_ini,
                      b_face_order,
                      b_face_n2o);

  cs_mesh_coherency_check();

  /* Field distribution
     ------------------ */

  // Vertex-centered fields.
  CS_FREE(vtx_n2o);

  // Cell-centered fields.
  // Note: the boundary coeffs will be exchanged after exchanging the boundary
  // faces.

  int n_fields = cs_field_n_fields();

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *field = cs_field_by_id(i);
    if (field->location_id != CS_MESH_LOCATION_CELLS) continue;

    // before resizing the field pointers, copy

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(field->location_id);

    cs_real_t **copy;
    CS_MALLOC(copy, field->n_time_vals, cs_real_t *);
    for (int j = 0; j < field->n_time_vals; j++) {
      CS_MALLOC(copy[j], n_cells_ini*field->dim, cs_real_t);
      memcpy(copy[j], field->vals[j], n_cells_ini*field->dim*sizeof(cs_real_t));

      field->_vals[j]->reshape(n_elts[2], field->dim);
      field->vals[j] = field->_vals[j]->data();
    }

    cs_real_t *recv_buf;
    CS_MALLOC(recv_buf, mesh->n_cells*field->dim, cs_real_t);

    for (int j = 0; j < field->n_time_vals; j++) {
      cs_all_to_all_copy_array(cd,
                               field->dim,
                               false,
                               copy[j],
                               recv_buf);

      CS_FREE(copy[j]);

      for (cs_lnum_t new_slot = 0; new_slot < mesh->n_cells; new_slot++) {
        cs_lnum_t old_slot = cell_n2o ? cell_n2o[new_slot] : new_slot;
        cs_lnum_t pre_slot = cell_order[old_slot];

        const cs_real_t *orig = recv_buf + field->dim*pre_slot;
        cs_real_t *dest = field->vals[j] + field->dim*new_slot;
        for (int dim = 0; dim < field->dim; dim++) {
          dest[dim] = orig[dim];
        }
      }
    }

    // Update val and val_pre pointer
    field->val = field->vals[0];
    if (field->n_time_vals > 1)
      field->val_pre = field->vals[1];

    CS_FREE(recv_buf);
    CS_FREE(copy);

    // exchange bc_coeffs
    cs_field_bc_coeffs_t *fc = field->bc_coeffs;

    if (!fc) continue;

    int i_mult, a_mult, b_mult;
    cs_field_get_bc_coeff_mult(field, &i_mult, &a_mult, &b_mult);

    cs_real_t *a_copy;
    cs_real_t *b_copy;
    CS_MALLOC(a_copy, n_b_faces_ini*a_mult, cs_real_t);
    CS_MALLOC(b_copy, n_b_faces_ini*b_mult, cs_real_t);

    int *i_copy;
    CS_MALLOC(i_copy, n_b_faces_ini*i_mult, int);
    memcpy(i_copy, fc->icodcl, n_b_faces_ini*i_mult*sizeof(int));
    CS_REALLOC(fc->icodcl, mesh->n_b_faces*i_mult, int);

    _exchange_and_order_b(bfd,
                          i_mult,
                          false,
                          i_copy,
                          fc->icodcl,
                          b_face_order,
                          b_face_n2o);

    CS_FREE(i_copy);

    // Non-interleaved bc codes
    cs_real_t *rcod;
    CS_MALLOC(rcod, n_b_faces_ini*a_mult, cs_real_t);

    memcpy(rcod, fc->rcodcl1, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(fc->rcodcl1, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *buf = fc->rcodcl1 + dim*mesh->n_b_faces;
      _exchange_and_order_b(bfd,
                            1,
                            false,
                            rcod + dim*n_b_faces_ini,
                            buf,
                            b_face_order,
                            b_face_n2o);
    }

    memcpy(rcod, fc->rcodcl2, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(fc->rcodcl2, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *buf = fc->rcodcl2 + dim*mesh->n_b_faces;
      _exchange_and_order_b(bfd,
                            1,
                            false,
                            rcod + dim*n_b_faces_ini,
                            buf,
                            b_face_order,
                            b_face_n2o);
    }

    memcpy(rcod, fc->rcodcl3, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(fc->rcodcl3, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *buf = fc->rcodcl3 + dim*mesh->n_b_faces;
      _exchange_and_order_b(bfd,
                            1,
                            false,
                            rcod + dim*n_b_faces_ini,
                            buf,
                            b_face_order,
                            b_face_n2o);
    }

    CS_FREE(rcod);

    _distribute_bc_coeff(bfd, n_b_faces_ini, a_mult, a_copy, &fc->a,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, b_mult, b_copy, &fc->b,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, a_mult, a_copy, &fc->af,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, b_mult, b_copy, &fc->bf,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, a_mult, a_copy, &fc->ad,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, b_mult, b_copy, &fc->bd,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, a_mult, a_copy, &fc->ac,
                         b_face_order, b_face_n2o);
    _distribute_bc_coeff(bfd, n_b_faces_ini, b_mult, b_copy, &fc->bc,
                         b_face_order, b_face_n2o);

    _distribute_bc_coeff(bfd, n_b_faces_ini, 1, a_copy, &fc->hint,
                         b_face_order, b_face_n2o);

    bool limited_val_f = fc->val_f == fc->val_f_lim;
    _distribute_bc_coeff(bfd, n_b_faces_ini, field->dim, a_copy, &fc->val_f,
                         b_face_order, b_face_n2o);
    if (limited_val_f) {
      fc->val_f_lim = fc->val_f;
    }
    else {
      _distribute_bc_coeff(bfd, n_b_faces_ini, field->dim, a_copy,
                           &fc->val_f_lim, b_face_order, b_face_n2o);
    }

    bool limited_flux = fc->flux == fc->flux_lim;
    _distribute_bc_coeff(bfd, n_b_faces_ini, field->dim, a_copy, &fc->flux,
                         b_face_order, b_face_n2o);
    if (limited_flux) {
      fc->flux_lim = fc->flux;
    }
    else {
      _distribute_bc_coeff(bfd, n_b_faces_ini, field->dim, a_copy,
                           &fc->flux_lim, b_face_order, b_face_n2o);
    }

    _distribute_bc_coeff(bfd, n_b_faces_ini, field->dim, a_copy,
                         &fc->val_f_pre, b_face_order, b_face_n2o);

    CS_FREE(a_copy);
    CS_FREE(b_copy);
  }

  cs_all_to_all_destroy(&cd);

  // Distribute the boundary fields here.

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *field = cs_field_by_id(i);
    if (field->location_id != CS_MESH_LOCATION_BOUNDARY_FACES) continue;

    cs_real_t **copy;
    CS_MALLOC(copy, field->n_time_vals*field->dim, cs_real_t *);
    for (int j = 0; j < field->n_time_vals; j++) {
      CS_MALLOC(copy[j], n_b_faces_ini*field->dim, cs_real_t);
      memcpy(copy[j], field->vals[j], n_b_faces_ini*field->dim*sizeof(cs_real_t));

      field->_vals[j]->reshape(mesh->n_b_faces, field->dim);
      field->vals[j] = field->_vals[j]->data();
    }

    cs_real_t *recv_buf;
    CS_MALLOC(recv_buf, mesh->n_b_faces*field->dim, cs_real_t);

    for (int j = 0; j < field->n_time_vals; j++) {
      cs_all_to_all_copy_array(bfd,
                               field->dim,
                               false,
                               copy[j],
                               recv_buf);

      CS_FREE(copy[j]);

      for (cs_lnum_t new_slot = 0; new_slot < mesh->n_b_faces; new_slot++) {
        cs_lnum_t old_slot = b_face_n2o ? b_face_n2o[new_slot] : new_slot;
        cs_lnum_t pre_slot = b_face_order[old_slot];

        const cs_real_t *orig = recv_buf + field->dim*pre_slot;
        cs_real_t *dest = field->vals[j] + field->dim*new_slot;

        for (int dim = 0; dim < field->dim; dim++) {
          dest[dim] = orig[dim];
        }
      }
    }

    field->val = field->vals[0];
    if (field->n_time_vals > 1)
      field->val_pre = field->vals[1];

    CS_FREE(recv_buf);
    CS_FREE(copy);
  }

  CS_FREE(b_face_order);
  cs_all_to_all_destroy(&bfd);

  // Distribute internal face-centered fields here

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *field = cs_field_by_id(i);
    if (field->location_id != CS_MESH_LOCATION_INTERIOR_FACES) continue;

    cs_real_t **copy;
    CS_MALLOC(copy, field->n_time_vals, cs_real_t *);
    for (int j = 0; j < field->n_time_vals; j++) {
      CS_MALLOC(copy[j], n_i_f_s*field->dim, cs_real_t);
      for (cs_lnum_t k = 0; k < n_i_f_s; k++) {
        cs_lnum_t f_id = i_face_list[k];
        for (int l = 0; l < field->dim; l++)
          copy[j][field->dim*k+l] = field->vals[j][field->dim*f_id+l];
      }

      field->_vals[j]->reshape(mesh->n_i_faces, field->dim);
      field->vals[j] = field->_vals[j]->data();
    }

    cs_real_t *recv_buf;
    CS_MALLOC(recv_buf, mesh->n_i_faces*field->dim, cs_real_t);

    for (int j = 0; j < field->n_time_vals; j++) {
      cs_all_to_all_copy_array(ifd,
                               field->dim,
                               false,
                               copy[j],
                               recv_buf);

      CS_FREE(copy[j]);

      for (cs_lnum_t new_slot = 0; new_slot < mesh->n_i_faces; new_slot++) {
        cs_lnum_t old_slot = i_face_n2o ? i_face_n2o[new_slot] : new_slot;
        cs_lnum_t pre_slot = i_face_order[old_slot];

        const cs_real_t *orig = recv_buf + field->dim*pre_slot;
        cs_real_t *dest = field->vals[j] + field->dim*new_slot;

        for (int dim = 0; dim < field->dim; dim++) {
          dest[dim] = orig[dim];
        }
      }
    }

    CS_FREE(recv_buf);
    CS_FREE(copy);

    field->val = field->vals[0];
    if (field->n_time_vals > 1)
      field->val_pre = field->vals[1];
  }

  CS_FREE(i_face_list);
  CS_FREE(i_face_order);
  cs_all_to_all_destroy(&ifd);

  // cs_preprocess_mesh_update_device();
  cs_gradient_free_quantities();
  cs_cell_to_vertex_free();
  cs_mesh_adjacencies_update_mesh();
  cs_matrix_update_mesh();

  CS_FREE(_dest_rank);
  CS_FREE(cell_order);
  CS_FREE(cell_n2o);
  CS_FREE(i_face_n2o);
  CS_FREE(b_face_n2o);
}

/*----------------------------------------------------------------------------*/
