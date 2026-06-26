/*============================================================================
 * Redistribution of mesh and field data.
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
#include "base/cs_order.h"
#include "base/cs_renumber_update.h"
#include "base/cs_field.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_sort.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "alge/cs_cell_to_vertex.h"
#include "alge/cs_gradient.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_redistribute.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compare cs_gnum_t (bsearch function).
 *
 * parameters:
 *   a <-- pointer to first number
 *   b <-- pointer to second number
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
 * Check for the existence of an array globally, then distribute it accordingly.
 *
 * For example, this is useful for boundary/AMR data that might not exist on
 * all procs.
 *
 * parameters:
 *   db          <-- pointer to the distributor.
 *   stride      <-- stride of the array to distribute.
 *   buf         <-> pointer to the array to redistribute.
 *----------------------------------------------------------------------------*/

template <typename T>
static void
_check_and_distribute_buffer(const cs_distributor_t   *db,
                             int                       stride,
                             T                        *buf[])
{
  assert(buf != nullptr);

  cs_lnum_t coeff_exists = (*buf != nullptr);
  cs_parall_max(1, CS_INT_TYPE, &coeff_exists);

  if (!coeff_exists)
    return;

  cs_distribute_buffer(db, stride, buf);
}

/*----------------------------------------------------------------------------
 * Distribute boundary conditions info and types.
 *
 * parameters:
 *   db          <-- pointer to the distributor.
 *----------------------------------------------------------------------------*/

static void
_distribute_bc_info_and_types(const cs_distributor_t   *db)
{
  cs_boundary_condition_pm_info_t *pm = cs_glob_bc_pm_info;

  cs_distribute_buffer(db, 1, &pm->izfppp);

  _check_and_distribute_buffer(db, 1, &pm->iautom);

  int **bc_type = cs_boundary_conditions_get_bc_type_addr();

  cs_distribute_buffer(db, 1, bc_type);

  cs_glob_bc_type = *bc_type;
}

/*----------------------------------------------------------------------------
 * Redistribute fields and bc_coeffs based on the mesh distributors.
 *
 * parameters:
 *   cd             <-- pointer to cell distributor.
 *   ifd            <-- pointer to interior face distributor.
 *   bfd            <-- pointer to boundary face distributor.
 *   bfd            <-- pointer to vertex distributor.
 *   n_cells_ini    <-- number of mesh cells before redistribution.
 *   n_i_faces_ini  <-- number of mesh interior faces before redistribution.
 *   n_b_faces_ini  <-- number of mesh boundary faces before redistribution.
 *   n_vertices     <-- number of mesh vertices before redistribution.
 *----------------------------------------------------------------------------*/

static void
_distribute_fields(const cs_distributor_t    *cd,
                   const cs_distributor_t    *ifd,
                   const cs_distributor_t    *bfd,
                   const cs_distributor_t    *vd,
                   const cs_lnum_t            n_cells_ini,
                   const cs_lnum_t            n_i_faces_ini,
                   const cs_lnum_t            n_b_faces_ini,
                   const cs_lnum_t            n_vertices_ini)
{
  cs_mesh_t *mesh = cs_glob_mesh;

  int n_fields = cs_field_n_fields();

  for (int i = 0; i < n_fields; i++) {

    cs_field_t *field = cs_field_by_id(i);

    const cs_distributor_t *db = nullptr;

    cs_mesh_location_type_t mlt = cs_mesh_location_get_type(field->location_id);
    cs_lnum_t n_elts_ini = 0;
    cs_lnum_t n_elts = 0;

    switch (mlt) {
      case CS_MESH_LOCATION_NONE:
        continue;

      case CS_MESH_LOCATION_CELLS:
        db = cd;
        n_elts_ini = n_cells_ini;
        n_elts = mesh->n_cells_with_ghosts;
        break;

      case CS_MESH_LOCATION_INTERIOR_FACES:
        db = ifd;
        n_elts_ini = n_i_faces_ini;
        n_elts = mesh->n_i_faces;
        break;

      case CS_MESH_LOCATION_BOUNDARY_FACES:
        db = bfd;
        n_elts_ini = n_b_faces_ini;
        n_elts = mesh->n_b_faces;
        break;

      case CS_MESH_LOCATION_VERTICES:
        db = vd;
        n_elts_ini = n_vertices_ini;
        n_elts = mesh->n_vertices;
        break;

      default:
        assert(0);
    }

    cs_assert(mlt == (cs_mesh_location_type_t)field->location_id);

    cs_real_t *copy = nullptr;
    CS_MALLOC(copy, n_elts_ini*field->dim, cs_real_t);

    for (int j = 0; j < field->n_time_vals; j++) {
      memcpy(copy, field->vals[j], n_elts_ini*field->dim*sizeof(cs_real_t));
      field->_vals[j]->reshape(n_elts, field->dim);
      field->vals[j] = field->_vals[j]->data();

      cs_distribute_buffer_allocated(db, field->dim, copy, field->vals[j]);
    }

    CS_FREE(copy);

    field->val = field->vals[0];
    if (field->n_time_vals > 1)
      field->val_pre = field->vals[1];

    cs_field_bc_coeffs_t *bcc = field->bc_coeffs;

    if (!bcc || mlt != CS_MESH_LOCATION_CELLS)
      continue;

    cs_lnum_t i_mult = 1, a_mult = 1, b_mult = 1;
    cs_field_get_bc_coeff_mult(field, &i_mult, &a_mult, &b_mult);

    cs_distribute_buffer(bfd, i_mult, &bcc->icodcl);

    // Note: the rcodcl family of buffers is non-interleaved.
    cs_real_t *rcod;
    CS_MALLOC(rcod, n_b_faces_ini*a_mult, cs_real_t);

    memcpy(rcod, bcc->rcodcl1, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(bcc->rcodcl1, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *dst = bcc->rcodcl1 + dim*mesh->n_b_faces;
      cs_distribute_buffer_allocated(bfd,
                                     1,
                                     rcod + dim*n_b_faces_ini,
                                     dst);
    }

    memcpy(rcod, bcc->rcodcl2, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(bcc->rcodcl2, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *dst = bcc->rcodcl2 + dim*mesh->n_b_faces;
      cs_distribute_buffer_allocated(bfd,
                                     1,
                                     rcod + dim*n_b_faces_ini,
                                     dst);
    }

    memcpy(rcod, bcc->rcodcl3, n_b_faces_ini*a_mult*sizeof(cs_real_t));
    CS_REALLOC(bcc->rcodcl3, mesh->n_b_faces*a_mult, cs_real_t);
    for (int dim = 0; dim < field->dim; dim++) {
      cs_real_t *dst = bcc->rcodcl3 + dim*mesh->n_b_faces;
      cs_distribute_buffer_allocated(bfd,
                                     1,
                                     rcod + dim*n_b_faces_ini,
                                     dst);
    }

    CS_FREE(rcod);

    _check_and_distribute_buffer(bfd, a_mult, &bcc->a);
    _check_and_distribute_buffer(bfd, b_mult, &bcc->b);

    _check_and_distribute_buffer(bfd, a_mult, &bcc->af);
    _check_and_distribute_buffer(bfd, b_mult, &bcc->bf);

    _check_and_distribute_buffer(bfd, a_mult, &bcc->ad);
    _check_and_distribute_buffer(bfd, b_mult, &bcc->bd);

    _check_and_distribute_buffer(bfd, a_mult, &bcc->ac);
    _check_and_distribute_buffer(bfd, b_mult, &bcc->bc);

    _check_and_distribute_buffer(bfd, 1, &bcc->h_int_tot);

    _check_and_distribute_buffer(bfd, field->dim, &bcc->val_f);
    _check_and_distribute_buffer(bfd, field->dim, &bcc->val_f_pre);

    _check_and_distribute_buffer(bfd, field->dim, &bcc->flux_diff);
  }

  _distribute_bc_info_and_types(bfd);
}

/*----------------------------------------------------------------------------
 * Free field gradients.
 *
 * Gradients are assumed to be computed on-the-fly when needed.
 *----------------------------------------------------------------------------*/

static void
_free_field_gradients(void)
{
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->grad)
      CS_FREE(f->grad);
  }
}

/*----------------------------------------------------------------------------
 * Generate random cell destination ranks (useful for stress testing).
 *
 * parameters:
 *   mesh          <-- pointer to the mesh
 *
 * returns:
 *   a random cell destination rank array.
 *----------------------------------------------------------------------------*/

static int *
_random_dest_rank(const cs_mesh_t  *mesh)
{
  MPI_Comm comm = cs_glob_mpi_comm;

  /* Generate and sync random seed. */

  static bool rng_on = false;
  if (!rng_on) {
    rng_on = true;

    unsigned int seed;

    if (cs_glob_rank_id <= 0) {
      seed = time(nullptr);
      //seed = 123456; // Fix the seed for debugging.
    }

    MPI_Bcast(&seed, 1, MPI_UINT32_T, 0, comm);

    srand(seed);
  }

  /* Allocate destination rank. */

  int *dest_rank = nullptr;
  CS_MALLOC(dest_rank, mesh->n_cells_with_ghosts, int);

  // Build the cell part-to-part distributor.

  const int n_max_tries = 50;
  cs_lnum_t g_min_n_cells = 0;

  for (int i = 0; i < n_max_tries && g_min_n_cells == 0; i++) {

    for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
      dest_rank[c_id] = rand() % cs_glob_n_ranks;

    cs_all_to_all_t *cd = cs_all_to_all_create(mesh->n_cells,
                                               0,
                                               nullptr,
                                               dest_rank,
                                               comm);

    g_min_n_cells = cs_all_to_all_n_elts_dest(cd);

    cs_all_to_all_destroy(&cd);

    cs_parall_min(1, CS_INT_TYPE, &g_min_n_cells);

  }

  cs_assert(g_min_n_cells > 0);

  cs_halo_sync_untyped(mesh->halo, CS_HALO_STANDARD, sizeof(int), dest_rank);

  return dest_rank;
}

/*----------------------------------------------------------------------------
 * Build a generic mesh elements distributor.
 *
 * This structure aims to distribute cells, interior/boundary faces and
 * vertices in a unified way, while minimizing duplicates in the send_list
 * array.
 *
 * The pointer to send_list can be set to null if all the elements are to be
 * distributed.
 *
 * After calling this function, if the pointer to send_list is not null,
 * the ownership of send_list is transferred to the distributor. Therefore, it
 * should not freed during the lifetime of the distributor.
 *
 * The gnum array corresponds to the elements' global ids. If send_list is not
 * null, the global id of the i-th element in send_list should correspond to
 * gnum[send_list[i]], and to gnum[i] otherwise.
 *
 * parameters:
 *   n_send              <-- number of elements to distribute
 *   send_list           <-- array of ids of the elements to send, of size
 *                           n_send, or nullptr.
 *   dest_ranks          <-- array of destination ranks of sent elements,
 *                           of size n_send.
 *   comm                <-- MPI communicator.
 *   gnum                <-- global ids the elements.
 *
 * returns:
 *   a part-to-part distributor.
 *----------------------------------------------------------------------------*/

static cs_distributor_t *
_distributor_create(cs_lnum_t           n_send,
                    cs_lnum_t          *send_list[],
                    const int           dest_ranks[],
                    const MPI_Comm      comm,
                    const cs_gnum_t     gnum[])
{
  cs_all_to_all_t *d = cs_all_to_all_create(n_send,
                                            0,
                                            nullptr,
                                            dest_ranks,
                                            comm);

  cs_lnum_t *_send_list = nullptr;

  if (send_list)
    _send_list = *send_list;

  cs_gnum_t *gnum_s = nullptr;
  CS_MALLOC(gnum_s, n_send, cs_gnum_t);

  for (cs_lnum_t i = 0; i < n_send; i++) {
    cs_lnum_t id = _send_list ? _send_list[i] : i;
    gnum_s[i] = gnum[id];
  }

  cs_gnum_t *gnum_r = cs_all_to_all_copy_array(d,
                                               1,
                                               false,
                                               gnum_s);

  CS_FREE(gnum_s);

  cs_lnum_t n_recv = cs_all_to_all_n_elts_dest(d);

  cs_lnum_t *recv_order = cs_order_gnum(nullptr, gnum_r, n_recv);

  cs_lnum_t *ordered_to_local = nullptr;
  CS_MALLOC(ordered_to_local, n_recv, cs_lnum_t);
  memset(ordered_to_local, -1, n_recv * sizeof(cs_lnum_t));

  cs_lnum_t n_uniq = 0;
  cs_gnum_t first = 0;

  for (cs_lnum_t i = 0; i < n_recv; i++) {
    cs_lnum_t ordered = recv_order[i];
    cs_gnum_t cur_gnum = gnum_r[ordered];
    assert(cur_gnum >= first);
    if (first != cur_gnum) {
      ordered_to_local[ordered] = n_uniq++;
      first = cur_gnum;
    }
  }

  CS_FREE(gnum_r);

  cs_distributor_t *db = nullptr;
  CS_MALLOC(db, 1, cs_distributor_t);

  db->d = d;
  db->n_send = n_send;
  db->n_recv = n_recv;
  db->n_uniq = n_uniq;
  db->send_list = _send_list;
  db->recv_order = recv_order;
  db->ordered_to_local = ordered_to_local;

  return db;
}

/*----------------------------------------------------------------------------
 * Create a vertex distributor.
 *
 * The vertices list and destination rank list are created based on the
 * cell-vertices connectivity and the cell destination ranks.
 *
 * The distributor guarantees the uniqueness of the vertices
 * post-redistribution.
 *
 * parameters:
 *   mesh              <-> pointer to the mesh
 *   cell_dest         <-- destination ranks of the cells.
 *   comm              <-- MPI communicator.
 *
 * returns:
 *   a vertex distributor.
 *----------------------------------------------------------------------------*/

static cs_distributor_t *
_create_vertex_distributor(const cs_mesh_t   *mesh,
                           const int          cell_dest[],
                           const MPI_Comm     comm)
{
  const cs_lnum_t n_vertices = mesh->n_vertices;

  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();
  assert(c2v->n_elts == mesh->n_cells);

  cs_adjacency_t *v2c = cs_adjacency_transpose(n_vertices, c2v);

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    for (cs_lnum_t i = v2c->idx[v_id]; i < v2c->idx[v_id+1]; i++) {
      v2c->ids[i] = cell_dest[v2c->ids[i]];
    }
  }

  cs_sort_indexed(n_vertices, v2c->idx, v2c->ids);

  cs_lnum_t n_elem = 0;

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    int prev_dest = -1;
    for (cs_lnum_t i = v2c->idx[v_id]; i < v2c->idx[v_id+1]; i++) {
      int dest = v2c->ids[i];
      assert(prev_dest <= dest);
      if (dest != prev_dest) {
        n_elem++;
        prev_dest = dest;
      }
    }
  }

  cs_lnum_t *v_list = nullptr;
  int *v_dests = nullptr;
  CS_MALLOC(v_list, n_elem, cs_lnum_t);
  CS_MALLOC(v_dests, n_elem, int);

  n_elem = 0;

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    int prev_dest = -1;
    for (cs_lnum_t i = v2c->idx[v_id]; i < v2c->idx[v_id+1]; i++) {
      int dest = v2c->ids[i];
      if (dest != prev_dest) {
        v_list[n_elem] = v_id;
        v_dests[n_elem] = dest;
        n_elem++;
        prev_dest = dest;
      }
    }
  }

  cs_adjacency_destroy(&v2c);

  cs_distributor_t *db = _distributor_create(n_elem,
                                             &v_list,
                                             v_dests,
                                             comm,
                                             mesh->global_vtx_num);

  cs_all_to_all_transfer_dest_rank(db->d, &v_dests);

  return db;
}

/*----------------------------------------------------------------------------
 * Create a boundary face distributor.
 *
 * parameters:
 *   mesh              <-> pointer to the mesh
 *   cell_dest         <-- destination ranks of the cells.
 *   comm              <-- MPI communicator.
 *
 * returns:
 *   a boundary face distributor.
 *----------------------------------------------------------------------------*/

static cs_distributor_t *
_create_b_face_distributor(const cs_mesh_t   *mesh,
                           const int          cell_dest_rank[],
                           const MPI_Comm     comm)
{
  int *dest_ranks = nullptr;
  CS_MALLOC(dest_ranks, mesh->n_b_faces, int);
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    dest_ranks[f_id] = cell_dest_rank[mesh->b_face_cells[f_id]];

  cs_distributor_t *db = _distributor_create(mesh->n_b_faces,
                                           nullptr,
                                           dest_ranks,
                                           comm,
                                           mesh->global_b_face_num);

  cs_all_to_all_transfer_dest_rank(db->d, &dest_ranks);

  return db;
}

/*----------------------------------------------------------------------------
 * Create an interior face distributor.
 *
 * The distributor guarantees the uniqueness of the interior faces
 * post-redistribution.
 *
 * parameters:
 *   mesh              <-> pointer to the mesh
 *   cell_dest         <-- destination ranks of the cells.
 *   comm              <-- MPI communicator.
 *
 * returns:
 *   an interior face distributor.
 *----------------------------------------------------------------------------*/

static cs_distributor_t *
_create_i_face_distributor(const cs_mesh_t   *mesh,
                           const int          cell_dest_rank[],
                           const MPI_Comm     comm)
{
  // We need the current ranks of the halo cells.

  int *cell_rank = nullptr;
  CS_MALLOC(cell_rank, mesh->n_cells_with_ghosts, int);
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
    cell_rank[c_id] = cs_glob_rank_id;
  cs_halo_sync_untyped(mesh->halo, CS_HALO_STANDARD, sizeof(int), cell_rank);

  cs_lnum_t *i_face_list = nullptr;
  int *i_face_dest = nullptr;

  CS_MALLOC(i_face_list, 2*mesh->n_i_faces, cs_lnum_t);
  CS_MALLOC(i_face_dest, 2*mesh->n_i_faces, int);

  cs_lnum_t n_send = 0;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {

    int prev_rank = -1;

    for (int i = 0; i < 2; i++) {
      int dest_rank = cell_dest_rank[mesh->i_face_cells[f_id][i]];
      if (dest_rank != prev_rank) {
        i_face_list[n_send] = f_id;
        i_face_dest[n_send] = dest_rank;
        n_send++;
        prev_rank = dest_rank;
      }
    }

  }

  CS_FREE(cell_rank);

  CS_REALLOC(i_face_list, n_send, cs_lnum_t);
  CS_REALLOC(i_face_dest, n_send, int);

  cs_distributor_t *ifd = _distributor_create(n_send,
                                            &i_face_list,
                                            i_face_dest,
                                            comm,
                                            mesh->global_i_face_num);

  cs_all_to_all_transfer_dest_rank(ifd->d, &i_face_dest);

  return ifd;
}

/*----------------------------------------------------------------------------
 * Create a cell distributor.
 *
 * The distributor guarantees the uniqueness of the interior faces
 * post-redistribution.
 *
 * parameters:
 *   mesh              <-> pointer to the mesh
 *   cell_dest         <-- destination ranks of the cells.
 *   comm              <-- MPI communicator.
 *
 * returns:
 *   a cell distributor.
 *----------------------------------------------------------------------------*/

static cs_distributor_t *
_create_cell_distributor(const cs_mesh_t   *mesh,
                         const int          cell_dest_rank[],
                         const MPI_Comm     comm)
{
  int *_cell_dest = nullptr;
  CS_MALLOC(_cell_dest, mesh->n_cells, int);
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
    _cell_dest[c_id] = cell_dest_rank[c_id];

  cs_distributor_t *cd = _distributor_create(mesh->n_cells,
                                             nullptr,
                                             _cell_dest,
                                             comm,
                                             mesh->global_cell_num);

  cs_all_to_all_transfer_dest_rank(cd->d, &_cell_dest);

  return cd;
}

#endif // defined(HAVE_MPI)

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a distributor structure.
 *
 * \param[in,out]  db       pointer of pointer to a cs_distributor_structure
 */
/*----------------------------------------------------------------------------*/

#if defined (HAVE_MPI)

void
cs_distributor_destroy(cs_distributor_t  **db)
{
  if (db == nullptr)
    return;

  cs_distributor_t *_db = *db;

  if (_db == nullptr)
    return;

  cs_all_to_all_destroy(&_db->d);
  CS_FREE(_db->send_list);
  CS_FREE(_db->recv_order);
  CS_FREE(_db->ordered_to_local);
  _db->n_send = -1;
  _db->n_recv = -1;
  _db->n_uniq = -1;

  CS_FREE(_db);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Redistribute mesh and field data based on a cell destination rank
 * map.
 *
 * If no cell map is given, a random one is created internally.
 *
 * The different distributors are returned if any of the pointers are
 * not null. In this case, it is up to the caller to free them whenever
 * they are no longer needed.
 *
 * \param[in]  cell_dest_rank  destination rank for each cell
 * \param[out] cell_db         pointer of pointer to cell redistributor
 * \param[out] i_face_db       pointer of pointer to interior face redistributor
 * \param[out] b_face_db       pointer of pointer to boundary face redistributor
 * \param[out] vertex_db       pointer of pointer to vertex redistributor
 */
/*----------------------------------------------------------------------------*/

void
cs_redistribute(const int           cell_dest_rank[],
                cs_distributor_t  **cell_db,
                cs_distributor_t  **i_face_db,
                cs_distributor_t  **b_face_db,
                cs_distributor_t  **vertex_db)
{
  if (cs_glob_n_ranks == 1)
    return;

  cs_mesh_t *mesh = cs_glob_mesh;

  /* Generate random destination rank if none is provided. */

  int *_dest_rank = nullptr;
  if (cell_dest_rank == nullptr) {
    _dest_rank = _random_dest_rank(mesh);
    cell_dest_rank = _dest_rank;
  }

  /* Free some elements which will need to be rebuilt when/if used.
     Halos are needed at some interior steps so destroyed later. */

  cs_mesh_free_rebuildable(mesh, false);
  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);

  _free_field_gradients();
  cs_gradient_free_quantities();
  cs_cell_to_vertex_free();

  MPI_Comm comm = cs_glob_mpi_comm;

  const cs_lnum_t n_cells_ini = mesh->n_cells;
  const cs_lnum_t n_b_faces_ini = mesh->n_b_faces;
  const cs_lnum_t n_i_faces_ini = mesh->n_i_faces;
  const cs_lnum_t n_vertices_ini = mesh->n_vertices;

  /* Create the vertex distributor before the mesh changes. */

  cs_distributor_t *vd = _create_vertex_distributor(mesh,
                                                    cell_dest_rank,
                                                    comm);

  /* Boundary faces */

  cs_distributor_t *bfd = _create_b_face_distributor(mesh,
                                                     cell_dest_rank,
                                                     comm);

  const cs_lnum_t n_b_faces = bfd->n_uniq;

  cs_distribute_buffer(bfd, 1, &mesh->global_b_face_num);

  _check_and_distribute_buffer(bfd, 1, &mesh->b_face_r_c_idx);

  cs_distribute_buffer(bfd, 1, &mesh->b_face_family);

  cs_gnum_t *g_b_face_cells = nullptr;
  CS_MALLOC(g_b_face_cells, n_b_faces_ini, cs_gnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces_ini; f_id++)
    g_b_face_cells[f_id] = mesh->global_cell_num[mesh->b_face_cells[f_id]];

  CS_FREE(mesh->b_face_cells);

  cs_distribute_buffer(bfd, 1, &g_b_face_cells);

  cs_gnum_t *g_b_face_vtx_lst = nullptr;
  CS_MALLOC(g_b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_gnum_t);
  for (cs_lnum_t i = 0; i < mesh->b_face_vtx_connect_size; i++)
    g_b_face_vtx_lst[i] = mesh->global_vtx_num[mesh->b_face_vtx_lst[i]];

  CS_FREE(mesh->b_face_vtx_lst);

  cs_distribute_buffer_indexed(bfd,
                               &mesh->b_face_vtx_idx,
                               &g_b_face_vtx_lst);

  mesh->b_face_vtx_connect_size = mesh->b_face_vtx_idx[n_b_faces];

  /* Interior faces */

  cs_distributor_t *ifd = _create_i_face_distributor(mesh,
                                                     cell_dest_rank,
                                                     comm);

  cs_lnum_t n_i_faces = ifd->n_uniq;

  _check_and_distribute_buffer(ifd, 1, &mesh->i_face_r_gen);

  cs_distribute_buffer(ifd, 1, &mesh->global_i_face_num);

  cs_distribute_buffer(ifd, 1, &mesh->i_face_family);

  cs_gnum_t *cell_gnum = cs_mesh_get_cell_gnum(mesh, 1);

  cs_gnum_t *g_i_face_cells = nullptr;
  CS_MALLOC(g_i_face_cells, 2*n_i_faces_ini, cs_gnum_t);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces_ini; f_id++) {
    for (int i = 0; i < 2; i++)
      g_i_face_cells[2*f_id+i] = cell_gnum[mesh->i_face_cells[f_id][i]];
  }

  CS_FREE(mesh->i_face_cells);
  CS_FREE(cell_gnum);

  cs_distribute_buffer(ifd, 2, &g_i_face_cells);

  cs_gnum_t *g_i_face_vtx_lst = nullptr;
  CS_MALLOC(g_i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_gnum_t);
  for (cs_lnum_t i = 0; i < mesh->i_face_vtx_connect_size; i++)
    g_i_face_vtx_lst[i] = mesh->global_vtx_num[mesh->i_face_vtx_lst[i]];

  CS_FREE(mesh->i_face_vtx_lst);

  cs_distribute_buffer_indexed(ifd,
                               &mesh->i_face_vtx_idx,
                               &g_i_face_vtx_lst);

  mesh->i_face_vtx_connect_size = mesh->i_face_vtx_idx[n_i_faces];

  /* Vertices */

  mesh->n_vertices = vd->n_uniq;

  cs_distribute_buffer(vd, 1, &mesh->global_vtx_num);

  cs_distribute_buffer(vd, 3, &mesh->vtx_coord);

  _check_and_distribute_buffer(vd, 1, &mesh->vtx_r_gen);

  /* Transform global to local face-vertices */

  // Interior faces.

  CS_MALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_lnum_t);

  for (cs_lnum_t i = 0; i < mesh->i_face_vtx_connect_size; i++) {
    cs_gnum_t gv = g_i_face_vtx_lst[i];
    void *found = bsearch(&gv, mesh->global_vtx_num, mesh->n_vertices,
                          sizeof(cs_gnum_t), _cmp_gnum);
    cs_assert(found);
    mesh->i_face_vtx_lst[i] = (cs_gnum_t *)found - mesh->global_vtx_num;
  }

  CS_FREE(g_i_face_vtx_lst);

  mesh->n_i_faces = n_i_faces;

  // Boundary faces.

  CS_MALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_lnum_t);

  for (cs_lnum_t i = 0; i < mesh->b_face_vtx_connect_size; i++) {
    cs_gnum_t gv = g_b_face_vtx_lst[i];
    void *found = bsearch(&gv, mesh->global_vtx_num, mesh->n_vertices,
                          sizeof(cs_gnum_t), _cmp_gnum);
    cs_assert(found);
    mesh->b_face_vtx_lst[i] = (cs_gnum_t *)found - mesh->global_vtx_num;
  }

  CS_FREE(g_b_face_vtx_lst);

  mesh->n_b_faces = n_b_faces;
  mesh->n_b_faces_all = mesh->n_b_faces;

  /* Cells */

  cs_distributor_t *cd = _create_cell_distributor(mesh,
                                                  cell_dest_rank,
                                                  comm);

  cs_distribute_buffer(cd, 1, &mesh->global_cell_num);

  cs_distribute_buffer(cd, 1, &mesh->cell_family);

  mesh->n_cells = cd->n_uniq;

  /* Apply global-to-local face-cells connectivity */

  // Boundary faces.

  CS_MALLOC(mesh->b_face_cells, mesh->n_b_faces, cs_lnum_t);

  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_gnum_t gc = g_b_face_cells[f_id];
    void *found = bsearch(&gc, mesh->global_cell_num, mesh->n_cells,
                          sizeof(cs_gnum_t), _cmp_gnum);
    cs_assert(found);
    mesh->b_face_cells[f_id] = (cs_gnum_t *)found - mesh->global_cell_num;
  }

  CS_FREE(g_b_face_cells);

  // Interior faces.

  CS_MALLOC(mesh->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int i = 0; i < 2; i++) {
      cs_gnum_t gc = g_i_face_cells[2*f_id+i];
      void *found = bsearch(&gc, mesh->global_cell_num, mesh->n_cells,
                            sizeof(cs_gnum_t), _cmp_gnum);
      mesh->i_face_cells[f_id][i] = found ?
                                    (cs_gnum_t *)found - mesh->global_cell_num :
                                    -1;
    }
  }

  CS_FREE(g_i_face_cells);

  /* Update the halo */

  cs_mesh_free_rebuildable(mesh, true);

  mesh->n_cells_with_ghosts = mesh->n_cells;

  cs_mesh_init_halo(mesh, nullptr, mesh->halo_type, -1, true);

  cs_mesh_update_partial();

  /* Distribute the fields. */

  _distribute_fields(cd,
                     ifd,
                     bfd,
                     vd,
                     n_cells_ini,
                     n_i_faces_ini,
                     n_b_faces_ini,
                     n_vertices_ini);

  if (i_face_db)
    *i_face_db = ifd;
  else
    cs_distributor_destroy(&ifd);

  if (b_face_db)
    *b_face_db = bfd;
  else
    cs_distributor_destroy(&bfd);

  if (vertex_db)
    *vertex_db = vd;
  else
    cs_distributor_destroy(&vd);

  if (cell_db)
    *cell_db = cd;
  else
    cs_distributor_destroy(&cd);

  CS_FREE(_dest_rank);
}

#endif // defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
