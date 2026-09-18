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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <cmath>
#include <chrono>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"

#include "alge/cs_gradient.h"
#include "alge/cs_matrix_default.h"

#include "base/cs_boundary_zone.h"
#include "base/cs_math.h"
#include "base/cs_post.h"
#include "base/cs_renumber.h"
#include "base/cs_volume_zone.h"

#include "mesh/cs_geom.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_algorithm.h"
#include "mesh/cs_mesh_boundary.h"
#include "mesh/cs_mesh_group.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_remove.h"
#include "mesh/cs_stl.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_cut.h"

/*----------------------------------------------------------------------------
 *  Global variables
 *----------------------------------------------------------------------------*/

cs_mesh_cut_options_t cs_glob_mesh_cut_options = {
  0.05,  /* poro_min */
  0.1    /* eps_corr_grad_lin */
};

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

constexpr cs_real_t _plane_tol = 1e-12;
constexpr int       _default_family_id = 1;
#define _DEBUG_ 0

/*============================================================================
 * Local structure definitions
 *============================================================================*/

enum {
  VTX_UNKNOWN = 0,
  VTX_FLUID = 1 << 0,
  VTX_SOLID = 1 << 1,
  VTX_INTERFACE = 1 << 2
};

/*------------------------------------------------------------------------------
 * This structure holds all the necessary buffers for the algorithm.
 * The buffers are allocated once, with enough memory to successfully handle
 * all the cells to be cut.
 * Unless specified, all the buffers are reset and filled on a cell-by-cell
 * basis.
 *
 * Members:
 *   sd       : signed distances of cell vertices to the cut plane.
 *   occurs   : number of occurrances of vertices in a polyline.
 *   e_stride : maximum number of edges per face.
 *   f_size   : maximum size of the face connectivity array induced by the
 *              cut.
 *   edges    : edge-vertex connectivity.
 *   e2f      : edge-face connectivity.
 *   xyz      : coordinates of the intersection points.
 *   indices  : indirection table, useful for eliminating duplicate vertices.
 *   compact  : maps unique vertices to their indices within the mesh.
 *   n_polys  : number of closing polygons created by the cutting algorithm.
 *              This variable gets incremented after every cell cut.
 *   polys    : the indices of the closing polygons are appended to this array.
 *----------------------------------------------------------------------------*/

struct _cut_data {
  cs_real_t    *sd;
  uint8_t      *occurs;
  cs_lnum_t     e_stride;
  cs_lnum_t     f_size;
  cs_lnum_2_t  *edges;
  cs_lnum_t    *e2f;
  cs_lnum_t    *f2e;
  cs_real_t    *xyz;
  cs_lnum_t    *indices;
  cs_lnum_t    *compact;
  cs_lnum_t     n_polys;
  cs_lnum_t    *polys;
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Remove invalid cells based on geometric criteria:
 *        cells in solid regions, cells with a low parent-to-child
 *        volume ratio, and highly warped cells.
 *
 * \param[in, out]  mesh            pointer to mesh structure
 * \param[in]       mq_old          pointer to primal mesh_quantities structure
 * \param[in]       n_cells_origin  previous number of cells
 * \param[in]       v2v             vertex to vertex connectivity
 * \param[in]       c2e             cell to edge connectivity
 * \param[in]       edge_v0         for each edge, index of vertex
 *                                  with smaller id
 * \param[in]       cell_flag       flag array marking cells disabled
 *                                  due to inconsistent vertex signs
 * \param[in]       n2o_cells       parent cell id for each cut cell
 * \param[in, out]  vertex_sign     sign of each vertex
 */
/*----------------------------------------------------------------------------*/

static void _remove_invalid_cells
  (cs_mesh_t                   *mesh,
   const cs_mesh_quantities_t  *mq_old,
   const cs_lnum_t              n_cells_origin,
   const cs_adjacency_t        *v2v,
   const cs_adjacency_t        *c2e,
   const cs_lnum_t              edge_v0[],
   const int                    cell_flag_inconsistency[],
   const cs_lnum_t              n2o_cells[],
   int                          vertex_sign[])
{
  const cs_lnum_t *v2v_ids = v2v->ids;
  const cs_lnum_t *c2e_idx = c2e->idx;
  const cs_lnum_t *c2e_ids = c2e->ids;

  /* Disabled cells nodes become IBM interface */

  for (cs_lnum_t c_id = 0; c_id < n_cells_origin; c_id++) {
    const cs_lnum_t s_id = c2e_idx[c_id];
    const cs_lnum_t e_id = c2e_idx[c_id+1];

    if (cell_flag_inconsistency[c_id] == 1) {

      for (cs_lnum_t _eidx = s_id; _eidx < e_id; _eidx++) {
        cs_lnum_t _edge_id = c2e_ids[_eidx];
        cs_lnum_t _v0 = edge_v0[_edge_id];
        cs_lnum_t _v1 = v2v_ids[_edge_id];

        vertex_sign[_v0] = VTX_INTERFACE;
        vertex_sign[_v1] = VTX_INTERFACE;
      }
    }
  }

  if (mesh->vtx_range_set == nullptr)
    mesh->vtx_range_set = cs_range_set_create(mesh->vtx_interfaces,
                                              nullptr,
                                              mesh->n_vertices,
                                              true, //false, // balance
                                              2,  // tr_ignore
                                              0); // g_id_base

  const cs_range_set_t *vrs = mesh->vtx_range_set;
  cs_datatype_t datatype = cs_datatype_from_type<int>();

  {
    cs_adjacency_t *v2v_new = cs_mesh_adjacency_v2v(mesh);
    const cs_lnum_t *v2v_new_idx = v2v_new->idx;
    const cs_lnum_t *v2v_new_ids = v2v_new->ids;

    cs_lnum_t sign_changed = 1;
    cs_lnum_t sign_changed_prev = -1;
    cs_lnum_t n_iter = 0;
    const cs_lnum_t n_iter_max = 1000;

    while (sign_changed > 0) {

      cs_range_set_zero_out_of_range(vrs, datatype, 1, vertex_sign);
      cs_range_set_sync(vrs, datatype, 1, vertex_sign);

      n_iter++;

      sign_changed = 0;
      for (cs_lnum_t v0 = 0; v0 < mesh->n_vertices; v0++) {

        for (cs_lnum_t i = v2v_new_idx[v0]; i < v2v_new_idx[v0+1]; i++) {
          cs_lnum_t v1 = v2v_new_ids[i];
          int s0 = vertex_sign[v0];
          int s1 = vertex_sign[v1];

          /* Propagates only fluid and solid sign */
          if (s0 == VTX_INTERFACE || s1 == VTX_INTERFACE)
            continue;

          if (s0 != VTX_UNKNOWN && s1 == VTX_UNKNOWN) {
            vertex_sign[v1] = s0;
            sign_changed += 1;
          }
          else if (s1 != VTX_UNKNOWN && s0 == VTX_UNKNOWN) {
            vertex_sign[v0] = s1;
            sign_changed += 1;
          }
        }
      }

      cs_parall_sum(1, CS_LNUM_TYPE, &sign_changed);

      if (n_iter >= n_iter_max && sign_changed == sign_changed_prev) {
        /* Detect stagnation only after a large number of iterations */
        bft_error(__FILE__, __LINE__, 0,
                  _("Vertex-sign propagation stopped after %d iterations: "
                    "no progress detected during the last iteration.\n"),
                    n_iter);
      }

      sign_changed_prev = sign_changed;

#if _DEBUG_
      if (cs_glob_rank_id == 0 || cs_glob_rank_id == -1)
        bft_printf("Propagation of sign: iteration %d : sign_changed = %d\n",
                   n_iter, sign_changed);
#endif

    }

    cs_adjacency_destroy(&v2v_new);

  }

  cs_range_set_sync(vrs, datatype, 1, vertex_sign);

  char *cell_flag;
  CS_MALLOC(cell_flag, mesh->n_cells, char);

  /* all is remove by default */
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
    cell_flag[c_id] = 1;

  for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[i][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[i][1];

    cs_lnum_t s_id = mesh->i_face_vtx_idx[i];
    cs_lnum_t e_id = mesh->i_face_vtx_idx[i+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t vtx_id = mesh->i_face_vtx_lst[j];
      int v_sign = vertex_sign[vtx_id];

      assert(   v_sign == VTX_FLUID
             || v_sign == VTX_SOLID
             || v_sign == VTX_INTERFACE);

      if (v_sign == VTX_FLUID) {
        if (c_id_0 < mesh->n_cells)
          cell_flag[c_id_0] = 0;
        if (c_id_1 < mesh->n_cells)
          cell_flag[c_id_1] = 0;
      }
    }
  }

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    cs_lnum_t c_id = mesh->b_face_cells[i];
    cs_lnum_t s_id = mesh->b_face_vtx_idx[i];
    cs_lnum_t e_id = mesh->b_face_vtx_idx[i+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t vtx_id = mesh->b_face_vtx_lst[j];
      int v_sign = vertex_sign[vtx_id];

      assert(   v_sign == VTX_FLUID
             || v_sign == VTX_SOLID
             || v_sign == VTX_INTERFACE);

      if (v_sign == VTX_FLUID) {
        cell_flag[c_id] = 0;
      }
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells_origin; c_id++) {
    if (cell_flag[c_id] == 0 && cell_flag_inconsistency[c_id] == 1)
      cell_flag[c_id] = 1;
  }

  cs_mesh_quantities_t *mq = cs_mesh_quantities_create();
  cs_mesh_quantities_compute(mesh, mq);

  cs_real_33_t *corr_grad_lin_inv = mq->corr_grad_lin;
  cs_real_t poro_min = cs_glob_mesh_cut_options.poro_min;
  double eps_corr_grad_lin = cs_glob_mesh_cut_options.eps_corr_grad_lin;
  cs_lnum_t n_bad_grad_lin = 0;

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    cs_lnum_t c_id_parent = n2o_cells[c_id];
    cs_real_t poro = mq->cell_vol[c_id] / mq_old->cell_vol[c_id_parent];

    if (poro < poro_min) {
      cell_flag[c_id] = 1;
    }

    double error_grad_lin = 0.;
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        const cs_real_t expected = (i == j) ? 1. : 0.;
        const cs_real_t diff = corr_grad_lin_inv[c_id][i][j] - expected;
        error_grad_lin += diff * diff;
      }
    }
    error_grad_lin = sqrt(error_grad_lin / 3.);

    if (error_grad_lin > eps_corr_grad_lin) {
      cell_flag[c_id] = 1;
      n_bad_grad_lin += 1;
#if _DEBUG_
      bft_printf("c_id = %d, err = %f\n"
                 "corr_grad_lin = %f %f %f\n"
                 "                %f %f %f\n"
                 "                %f %f %f\n\n",
                 c_id, error_grad_lin,
                 corr_grad_lin_inv[c_id][0][0], corr_grad_lin_inv[c_id][0][1],
                 corr_grad_lin_inv[c_id][0][2], corr_grad_lin_inv[c_id][1][0],
                 corr_grad_lin_inv[c_id][1][1], corr_grad_lin_inv[c_id][1][2],
                 corr_grad_lin_inv[c_id][2][0], corr_grad_lin_inv[c_id][2][1],
                 corr_grad_lin_inv[c_id][2][2]);
#endif
    }
  }

  bft_printf(" %d cells are removed by the corr_grad_lin criteria\n",
             n_bad_grad_lin);

  cs_mesh_remove_cells(mesh, cell_flag, "auto:closing_polygons");
  /* Mark for re-partitioning */
  mesh->modified |= CS_MESH_MODIFIED_BALANCE;

  cs_mesh_quantities_destroy(mq);

  /* Remove isolated cells */

  cs_mesh_adjacencies_initialize();
  cs_mesh_adjacencies_update_mesh();
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_lnum_t *c2c_idx = ma->cell_cells_idx;

  cs_lnum_t n_isolated_cells = 0;

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    const cs_lnum_t n_internal_faces = e_id_i - s_id_i;

    if (n_internal_faces > 0) {
      cell_flag[c_id] = 0;
    }
    else {
      cell_flag[c_id] = 1;
      n_isolated_cells++;
    }
  }

  bft_printf(" %d isolated cells are removed\n",
             n_isolated_cells);

  cs_mesh_remove_cells(mesh, cell_flag, "auto:closing_polygons");

  cs_mesh_adjacencies_finalize();
  CS_FREE(cell_flag);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Update a global numbering array.
 *
 * Since the algorithm only creates new unique entities, a simple scan is
 * enough.
 *
 * \param[in, out]  global_num  the global numbering array to update.
 * \param[in]       n_local     the local number of entities post-cut.
 * \param[in]       delta       the number of new entities created by the cut.
 * \param[in, out]  n_global    the global number of entities post-cut.
 */
/*----------------------------------------------------------------------------*/

static void
_update_global_num(cs_gnum_t *global_num[],
                   cs_lnum_t  n_local,
                   cs_lnum_t  delta,
                   cs_gnum_t *n_global)
{
  cs_gnum_t local_new = (cs_gnum_t)delta;
  cs_gnum_t scan_new = local_new;
  cs_gnum_t total_new = local_new;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    MPI_Scan(&local_new, &scan_new, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&local_new, &total_new, 1, CS_MPI_GNUM, MPI_SUM,
        cs_glob_mpi_comm);
  }

#endif // defined(HAVE_MPI)

  if (*global_num) {
    CS_REALLOC(*global_num, n_local, cs_gnum_t);
    cs_gnum_t *_global_num = *global_num;
    cs_gnum_t start = scan_new - local_new;

    for (cs_lnum_t i = 0; i < delta; i++)
      _global_num[n_local-delta+i] = *n_global + start + i + 1;
  }

  *n_global = *n_global + total_new;
}

static void
_update_global_vertices(cs_mesh_t *mesh, cs_lnum_t n_new_vertices)
{
  _update_global_num(&mesh->global_vtx_num,
                     mesh->n_vertices,
                     n_new_vertices,
                     &mesh->n_g_vertices);
}

static void
_update_global_b_faces(cs_mesh_t *mesh, cs_lnum_t n_new_b_faces)
{
  _update_global_num(&mesh->global_b_face_num,
                     mesh->n_b_faces,
                     n_new_b_faces,
                     &mesh->n_g_b_faces);
}

static void
_update_global_i_faces(cs_mesh_t *mesh, cs_lnum_t n_new_i_faces)
{
  _update_global_num(&mesh->global_i_face_num,
                     mesh->n_i_faces,
                     n_new_i_faces,
                     &mesh->n_g_i_faces);
}

static void
_update_global_cells(cs_mesh_t *mesh, cs_lnum_t n_new_cells)
{
  _update_global_num(&mesh->global_cell_num,
                     mesh->n_cells,
                     n_new_cells,
                     &mesh->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Update parallel mesh structures and counts.
 *
 * Since the algorithm only creates new unique entities, a simple scan is
 * enough.
 *
 * \param[in, out]   mesh                    pointer to mesh structure.
 * \param[in]        update_global_vertices  update global vertex information
 * \param[in]        n_new_vertices          new number of local vertices
 * \param[in]        update_global_b_faces   update global boundary face
 *                                           information
 * \param[in]        n_new_b_faces           new number of local boundary faces
 * \param[in]        update_global_i_faces   update global internal face
 *                                           information
 * \param[in]        n_new_i_faces           new number of local internal faces
 * \param[in]        n_new_cells             new number of local cells
 */
/*----------------------------------------------------------------------------*/

static void
_update_parallelism(cs_mesh_t *mesh,
                    bool       update_global_vertices,
                    cs_lnum_t  n_new_vertices,
                    bool       update_global_b_faces,
                    cs_lnum_t  n_new_b_faces,
                    bool       update_global_i_faces,
                    cs_lnum_t  n_new_i_faces,
                    cs_lnum_t  n_new_cells)
{
  mesh->n_cells_with_ghosts = mesh->n_cells;

  if (update_global_vertices)
    _update_global_vertices(mesh, n_new_vertices);

  if (update_global_b_faces)
    _update_global_b_faces(mesh, n_new_b_faces);

  if (update_global_i_faces)
    _update_global_i_faces(mesh, n_new_i_faces);

  _update_global_cells(mesh, n_new_cells);

  if (mesh->n_domains > 1 || mesh->n_init_perio > 0) {
    cs_halo_type_t halo_type = mesh->halo_type;
    cs_mesh_builder_t *mb = (mesh == cs_glob_mesh) ?
                            cs_glob_mesh_builder :
                            nullptr;
    cs_mesh_init_halo(mesh, mb, halo_type, -1, true);
  }
}

/*------------------------------------------------------------------------------
 * \brief Reverse the order of the elements in the range [arr, arr+n]
 *----------------------------------------------------------------------------*/

static void
_reverse_array(cs_lnum_t *arr, cs_lnum_t n)
{
  cs_lnum_t start = 0, end = n-1, temp;
  while (start < end) {
    temp = arr[start];
    arr[start] = arr[end];
    arr[end] = temp;
    start++;
    end--;
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Attempt to determine the correct orientation of a new cycle contour
 *        by comparing it to an existing sub-face contour (reference).
 *
 * This function checks whether the new contour (fv_new)
 * should be reversed to ensure consistent orientation with a reference
 * sub-face contour. It looks for at least one shared edge between
 * the two contours and decides the orientation based on the direction
 * of that common edge.
 *
 * Behavior:
 *   - If a common edge is found in the same circulation (p->q matches _p->_q):
 *     --> the new contour should be reversed
 *   - If a common edge is found in the opposite direction (p->q matches _q->_p):
 *     --> the new contour is not reversed
 *
 * The function modifies fv_new in place only if reversal is required.
 *
 * \param[in]     reorient       logical flag that inverts the reversal logic:
 *                               - false: same-circulation edge -> reverse the cycle
 *                               - true:  same-direction edge   -> keep the cycle
 * \param[in,out] fv_new         array of vertex ids forming the new cycle contour
 *                               (closed: fv_new[0] assumed == fv_new[ne_new-1]
 *                               or handled cyclically)
 * \param[in]     ne_new         number of edges (and vertices) in the new contour
 * \param[in]     sfv            reference sub-face vertex list (closed contour)
 * \param[in]     ne_sf          number of edges (and vertices) in the reference

 * \return                      true if the new contour need to be reversed,
 *                              false otherwise
 */
/*----------------------------------------------------------------------------*/

static bool inline
_reorient_cycle_contour(const bool       reorient,
                        cs_lnum_t       *fv_new,
                        const int        ne_new,
                        const cs_lnum_t *sfv,
                        const int        ne_sf)
{
  bool reverse = false;
  bool found = false;

  for (int j = 0; j < ne_sf; ++j) {
    cs_lnum_t _p = sfv[j];
    cs_lnum_t _q = sfv[(j+1)%ne_sf];

    for (int k = 0; k < ne_new; ++k) {
      cs_lnum_t p = fv_new[k];
      cs_lnum_t q = fv_new[(k+1)%ne_new];

      if (p == _p && q == _q) {
        found = true;

        if (!(reorient)) {
          reverse = true;
        }
        else {
          reverse = false;
        }
        break;
      }
      if (p == _q && q == _p) {
        found = true;

        if (!(reorient)) {
          reverse = false;
        }
        else {
          reverse = true;
        }
        break;
      }
    }
    if (found)
      break;
  }

  if (reverse) {
    _reverse_array(fv_new+1, ne_new-1);
  }

  return found;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Create new cells and immersed faces by closing intersection cycles
 *        inside cut cells.
 *
 * For each selected cut cell:
 *   - Reorders light edges into one or more closed cycles (polygonal
 *      contours) that represent the intersection polygons between the cutting
 *      plane(s) and the cell.
 *   - For each detected cycle:
 *      - Creates a new immersed interior face using the intersection points
 *        along the cycle
 *      - Assigns connectivity and cell adjacency for the new immersed face
 *   - Updates the cell adjacency of existing sub-faces (both interior and
 *     boundary) so that sub-faces lying inside the new (positive-side) cell
 *     point to the newly created cell instead of the original one.
 *
 *   Light edges represent the contour visible from the positive side.
 *   Only "light" edges are used to build cycles
 *   One new cell is created per cycle (can lead to multiple new
 *                                      cells per original cell)
 *
 * \param[in,out] mesh                      pointer to the mesh structure
 * \param[in]     n_i_faces_old             number of interior faces before
 *                                          any cutting
 * \param[in,out] n_new_cells               counter: number of newly created cells
 * \param[in]     n_cut_cells               number of cells selected for cutting
 * \param[in]     sel_cells                 list of cell ids to be processed
 * \param[in]     c2f                       cell to face adjacency
 * \param[in]     vertex_sign               vertex sign
 * \param[in]     e_v_idx                   start global id of new vertices
 *                                          created on each edge
 * \param[in]     cell_flag                 flag array marking cells disabled
 *                                          due to inconsistent vertex signs
 * \param[in]     b_face_o2n_idx            sub-face cumulative count per
 *                                          original boundary face
 * \param[in]     i_face_o2n_idx            sub-face cumulative count per
 *                                          original interior face
 * \param[in]     b_sub_face_vtx_idx        vertex index for boundary sub-faces
 * \param[in]     b_sub_face_vtx_lst        vertex list for boundary sub-faces
 * \param[in]     i_sub_face_vtx_idx        vertex index for interior sub-faces
 * \param[in]     i_sub_face_vtx_lst        vertex list for interior sub-faces
 * \param[in]     light_edge_in_b_face_idx  start index of light edges per
 *                                          boundary face
 * \param[in]     light_edge_in_i_face_idx  start index of light edges per
 *                                          interior face
 * \param[in]     light_edge_in_b_face      [edge_a, edge_b] pairs for light
 *                                           edges on boundary faces
 * \param[in]     light_edge_in_i_face      [edge_a, edge_b] pairs for light
 *                                           edges on interior faces
 * \param[in,out] n_i_faces_tot             total number of interior faces
 *                                          after adding immersed ones
 * \param[in,out] n_i_face_vtx_tot          total number of vertex entries in
 *                                          all interior faces
 * \param[out]    _b_face_cells             updated cell adjacency for all
 *                                          boundary sub-faces
 * \param[out]    _i_face_cells             updated cell adjacency pairs for all
 *                                          interior sub-faces
 */
/*----------------------------------------------------------------------------*/

static void
_cut_cells(cs_mesh_t             *mesh,
           const cs_lnum_t        n_i_faces_old,
           cs_lnum_t             &n_new_cells,
           const cs_lnum_t        n_cut_cells,
           const cs_lnum_t       *sel_cells,
           const cs_adjacency_t  *c2f,
           const int              vertex_sign[],
           const cs_lnum_t        e_v_idx[],
           const int              cell_flag[],
           const cs_lnum_t        b_face_o2n_idx[],
           const cs_lnum_t        i_face_o2n_idx[],
           const cs_lnum_t        b_sub_face_vtx_idx[],
           const cs_lnum_t        b_sub_face_vtx_lst[],
           const cs_lnum_t        i_sub_face_vtx_idx[],
           const cs_lnum_t        i_sub_face_vtx_lst[],
           const cs_lnum_t        light_edge_in_b_face_idx[],
           const cs_lnum_t        light_edge_in_i_face_idx[],
           const cs_lnum_2_t      light_edge_in_b_face[],
           const cs_lnum_2_t      light_edge_in_i_face[],
           cs_lnum_t              n2o_cells[],
           cs_lnum_t             &n_i_faces_tot,
           cs_lnum_t             &n_i_face_vtx_tot,
           cs_lnum_t              _b_face_cells[],
           cs_lnum_t              _i_face_cells[][2])
{
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_lnum_t *c2f_idx = c2f->idx;
  const cs_lnum_t *c2f_ids = c2f->ids;

  for (cs_lnum_t c_id_loc = 0; c_id_loc < n_cut_cells; c_id_loc++) {
    cs_lnum_t c_id = sel_cells[c_id_loc];

    if (cell_flag[c_id] == 1)
      continue;

    int n_edge_in_cell = 0;

    const cs_lnum_t s_id_c = c2f_idx[c_id];
    const cs_lnum_t e_id_c = c2f_idx[c_id+1];
    const cs_lnum_t n_faces = e_id_c - s_id_c;

    cs_lnum_t max_size_select = n_faces*16;
    cs_lnum_2_t *edge_in_cell;
    CS_MALLOC(edge_in_cell, max_size_select, cs_lnum_2_t);

    /* Collect light edges */

    for (cs_lnum_t cidx = s_id_c; cidx < e_id_c; cidx++) {
      const cs_lnum_t g_face_id = c2f_ids[cidx];
      cs_lnum_t n_edge_in_face = 0;

      const bool is_internal = (g_face_id < n_i_faces_old);

      const cs_lnum_t face_id = is_internal ?
                                g_face_id :
                                g_face_id - n_i_faces_old;

      const cs_lnum_t *light_edge_in_face_idx
        = (is_internal) ? light_edge_in_i_face_idx : light_edge_in_b_face_idx;

      const cs_lnum_2_t *light_edge_in_face = (is_internal) ?
                                               light_edge_in_i_face:
                                               light_edge_in_b_face;

      cs_lnum_t le0 = light_edge_in_face_idx[face_id];
      cs_lnum_t le1 = light_edge_in_face_idx[face_id+1];

      for (cs_lnum_t i = le0; i < le1; i++) {
        cs_lnum_t v0 = light_edge_in_face[i][0];
        cs_lnum_t v1 = light_edge_in_face[i][1];

        assert(n_edge_in_cell + n_edge_in_face < max_size_select);

        edge_in_cell[n_edge_in_cell + n_edge_in_face][0] = v0;
        edge_in_cell[n_edge_in_cell + n_edge_in_face][1] = v1;

        n_edge_in_face++;

      }

      n_edge_in_cell += n_edge_in_face;
    }

    /* Reorder the light edges to construct disjoint closed cycles */

    // Reference edge to start is the first */

    cs_lnum_t *reorder_edge_id;
    CS_MALLOC(reorder_edge_id, n_edge_in_cell, cs_lnum_t);

    cs_lnum_t edges_used_tot = 0;
    cs_lnum_t out_pos = 0;
    cs_lnum_t n_cycle = 0;

    bool *used;
    CS_MALLOC(used, n_edge_in_cell, bool);
    for (cs_lnum_t j = 0; j < n_edge_in_cell; j++) {
      used[j] = false;
    }

    /* Start index of each cycle in reorder_edge_id */
    const int max_n_cycle = n_edge_in_cell;
    int *cycle_idx;
    CS_MALLOC(cycle_idx, max_n_cycle, int);

    /* Cycle_length */
    int *cycle_len;
    CS_MALLOC(cycle_len, n_edge_in_cell, int);

    while (edges_used_tot < n_edge_in_cell) {

      // Find j0 the reference edge_id not used for the new cycle
      cs_lnum_t j0 = -1;
      for (cs_lnum_t j = 0; j < n_edge_in_cell; j++) {
        if (!used[j]) {
          j0 = j;
          break;
        }
      }
      if (j0 == -1) break; // all is used --> finish

      cycle_idx[n_cycle] = out_pos;

      // Initialize the cycle with this reference edge */
      cs_lnum_t start = edge_in_cell[j0][0];
      cs_lnum_t next  = edge_in_cell[j0][1];
      reorder_edge_id[out_pos] = start;
      reorder_edge_id[out_pos+1] = next;
      out_pos += 2;
      assert(out_pos <= n_edge_in_cell);
      cycle_len[n_cycle] = 2;

      used[j0] = true;
      edges_used_tot++;

      cs_lnum_t last = next;
      bool closed = false;

      // Propagates until closing (back to start) */

      while (!closed) {
        bool progressed = false;

        for (cs_lnum_t j = 0; j < n_edge_in_cell; j++) {
          if (used[j]) continue;

          cs_lnum_t e_id0 = edge_in_cell[j][0];
          cs_lnum_t e_id1 = edge_in_cell[j][1];

          if (!(e_id0 == last || e_id1 == last))
            continue;

          cs_lnum_t e_id_next = (e_id0 == last) ? e_id1 : e_id0;

          /* edge used */
          used[j] = true;
          edges_used_tot++;

          if (e_id_next == start) {
            /* cycle is closed */
            closed = true;
          }
          else {
            /* Add the new vertex in the current cycle */
            reorder_edge_id[out_pos] = e_id_next;
            out_pos++;
            assert(out_pos <= n_edge_in_cell);
            cycle_len[n_cycle]++;
            last = e_id_next;
          }

          progressed = true;
          break;

        }

        if (!progressed) {
          // Edge missing ?
          bft_error(__FILE__, __LINE__, 0,
                    _("Not success to closed the cycle"));
        }

      } /* Cycle is closed */

      n_cycle++;
    }

#if _DEBUG_
    bft_printf("c_id = %d, n_cycle = %d\n", c_id, n_cycle);
#endif

    /* Check if sum of cycle_len = n_edge_in_cell */
    int sum = 0;
    for (cs_lnum_t cycle_id = 0; cycle_id < n_cycle; cycle_id++) {
      sum += cycle_len[cycle_id];
    }
    assert(sum == n_edge_in_cell);

    for (cs_lnum_t cycle_id = 0; cycle_id < n_cycle; cycle_id++) {
      cs_lnum_t start_vtx = cycle_idx[cycle_id];
      cs_lnum_t cycle_length = cycle_len[cycle_id];

      /* One new cell for each cycle */
      cs_lnum_t c_id_new = n_cells + n_new_cells;

      int vtx_ids[20];
      assert(cycle_length < 20);

      for (cs_lnum_t j = 0; j < cycle_length; j++) {
        cs_lnum_t edge_id = reorder_edge_id[start_vtx + j];

#if _DEBUG_
        bft_printf("cycle_id = %d, j =%d, reorder_edge_id = %d, "
                   "vtx_id = %d, c_id_new = %d\n",
                   cycle_id, j, reorder_edge_id[start_vtx + j],
                   e_v_idx[reorder_edge_id[start_vtx + j]], c_id_new);
#endif

        const cs_lnum_t v_id = e_v_idx[edge_id];
        vtx_ids[j] = e_v_idx[edge_id];

        mesh->i_face_vtx_lst[n_i_face_vtx_tot + j] = v_id;
      }

      /* Update i_face_cell and b_face_cell
         Each cycle --> 1 new cell */

      for (cs_lnum_t cidx = s_id_c; cidx < e_id_c; cidx++) {
        const cs_lnum_t g_f_id = c2f_ids[cidx];

        const bool is_internal = (g_f_id < n_i_faces_old);

        const cs_lnum_t f_id = is_internal ?
                               g_f_id : g_f_id - n_i_faces_old;

        /* Reorient if normal towards to interior */
        bool reorient = false;
        if (is_internal)
          reorient = (i_face_cells[f_id][0] == c_id) ? false : true;

        const cs_lnum_t *face_o2n_idx = (is_internal) ?
                                         i_face_o2n_idx:
                                         b_face_o2n_idx;

        cs_lnum_t s_id = face_o2n_idx[f_id];
        cs_lnum_t e_id = face_o2n_idx[f_id+1];

        /* loop on sub faces */

        for (cs_lnum_t sub_f_id = s_id; sub_f_id < e_id; sub_f_id++) {

          const cs_lnum_t *sub_face_vtx_idx = is_internal ?
                                              i_sub_face_vtx_idx :
                                              b_sub_face_vtx_idx;

          const cs_lnum_t *sub_face_vtx_lst = is_internal ?
                                              i_sub_face_vtx_lst :
                                              b_sub_face_vtx_lst;

          const cs_lnum_t v0 = sub_face_vtx_idx[sub_f_id];
          const cs_lnum_t v1 = sub_face_vtx_idx[sub_f_id + 1];

          bool sub_face_in_cycle = false;
          bool process_reverse = false;
          bool all_vtx_sign_is_zero_or_one = true;
          bool all_vtx_sign_is_one = true;
          bool all_vtx_sign_is_minus_one = true;
          bool all_vtx_sign_is_zero_or_minus_one = true;
          bool all_vtx_sign_is_zero = true;

          for (cs_lnum_t k = v0; k < v1; k++) {
            const cs_lnum_t sv = sub_face_vtx_lst[k];
            const int sign = vertex_sign[sv];

            if (!(sign == VTX_INTERFACE || sign == VTX_SOLID)) {
              all_vtx_sign_is_zero_or_one = false;
            }
            if (!(sign == VTX_FLUID)) {
              all_vtx_sign_is_minus_one = false;
            }
            if (!(sign == VTX_SOLID)) {
              all_vtx_sign_is_one = false;
            }
            if (!(sign == VTX_INTERFACE || sign == VTX_FLUID)) {
              all_vtx_sign_is_zero_or_minus_one = false;
            }
            if (!(sign == VTX_INTERFACE)) {
              all_vtx_sign_is_zero = false;
            }

            if (!sub_face_in_cycle) {
              for (cs_lnum_t j = 0; j < cycle_length; j++) {
                if (vtx_ids[j] == sv) {
                  sub_face_in_cycle = true;
                  break;
                }
              }
            }

          }

#if _DEBUG_
          bft_printf("cycle = %d, c_id = %d, face_id = %d, sub_f_id = %d, "
                     "all_vtx_sign_is_zero_or_one = %d, one = %d, "
                     "minus one = %d, sub_face_in_cycle = %d \n",
                     cycle_id, c_id, g_f_id, sub_f_id,
                     all_vtx_sign_is_zero_or_one,
                     all_vtx_sign_is_one,
                     all_vtx_sign_is_minus_one,
                     sub_face_in_cycle);
          bft_printf("sub_face_vtx_lst --> ");
          for (cs_lnum_t k = v0; k < v1; k++)
            bft_printf("%d ", sub_face_vtx_lst[k]);
          bftprintf("\n");
#endif

          /* Compare the circulation of the cycle with the sub face.
             For closed contour, the circulation need to be opposite */
          if (   sub_face_in_cycle
              && !(process_reverse)
              && all_vtx_sign_is_zero_or_minus_one) {
              process_reverse
                = _reorient_cycle_contour(reorient,
                                          mesh->i_face_vtx_lst+n_i_face_vtx_tot,
                                          cycle_length,
                                          sub_face_vtx_lst + v0,
                                          v1-v0);
          }

          if (   sub_face_in_cycle
              && all_vtx_sign_is_zero_or_one
              && !(all_vtx_sign_is_zero)) {
#if _DEBUG_
            bft_printf("--> Update new cell :\n");
#endif
            if (is_internal) {
              if (i_face_cells[f_id][0] == c_id)
                _i_face_cells[sub_f_id][0] = n_cells + n_new_cells;
              else if (i_face_cells[f_id][1] == c_id)
                _i_face_cells[sub_f_id][1] = n_cells + n_new_cells;
#if _DEBUG_
              bft_printf("ifacecell = %d %d, _i_face_cells = %d %d\n",
                         i_face_cells[f_id][0], i_face_cells[f_id][1],
                         _i_face_cells[sub_f_id][0], _i_face_cells[sub_f_id][1]);
#endif
            }
            else {
              _b_face_cells[sub_f_id] = n_cells + n_new_cells;
#if _DEBUG_
              bft_printf("bfacecell = %d, _b_face_cells = %d\n",
                         mesh->b_face_cells[f_id], _b_face_cells[sub_f_id]);
#endif
            }
          }
          else {

            if (sub_face_in_cycle) {
#if _DEBUG_
              bft_printf("--> sub_face is in cycle, Keep standard cell :\n");
#endif
              if (is_internal) {
#if _DEBUG_
                bft_printf("f_id = %d, sub_f_id = %d, i_face_cells = %d %d,"
                           " _i_face_cells = %d %d\n",
                           f_id, sub_f_id,
                           i_face_cells[f_id][0], i_face_cells[f_id][1],
                           _i_face_cells[sub_f_id][0], _i_face_cells[sub_f_id][1]);
#endif
              }
              else {
#if _DEBUG_
                bft_printf("bfacecell = %d, _b_face_cells = %d\n",
                           mesh->b_face_cells[f_id], _b_face_cells[sub_f_id]);
#endif
              }
            }
            else {
              if (all_vtx_sign_is_one) {
#if _DEBUG_
                bft_printf("--> sub_face is NOT in cycle, "
                           "isolated face with vtx sign = 1\n");
#endif
                if (is_internal) {
                  if (i_face_cells[f_id][0] == c_id)
                    _i_face_cells[sub_f_id][0] = n_cells + n_new_cells;
                  else if (i_face_cells[f_id][1] == c_id)
                    _i_face_cells[sub_f_id][1] = n_cells + n_new_cells;
#if _DEBUG_
                  bft_printf("f_id = %d, sub_f_id = %d, "
                             "ifacecell = %d %d, _i_face_cells = %d %d\n",
                             f_id, sub_f_id,
                             i_face_cells[f_id][0], i_face_cells[f_id][1],
                             _i_face_cells[sub_f_id][0],
                             _i_face_cells[sub_f_id][1]);
#endif
                }
                else {
                  _b_face_cells[sub_f_id] = n_cells + n_new_cells;
#if _DEBUG_
                  bft_printf("bfacecell = %d, _b_face_cells = %d\n",
                             mesh->b_face_cells[f_id], _b_face_cells[sub_f_id]);
#endif

                }
              }
              else if (all_vtx_sign_is_minus_one) {
#if _DEBUG_
                bft_printf("--> sub_face is NOT in cycle, "
                       "isolated face with vtx sign = -1\n");
#endif
                if (is_internal) {
#if _DEBUG_
                  bft_printf("f_id = %d, sub_f_id = %d, "
                             "i_face_cells = %d %d, _i_face_cells = %d %d\n",
                             f_id, sub_f_id,
                             i_face_cells[f_id][0], i_face_cells[f_id][1],
                             _i_face_cells[sub_f_id][0],
                             _i_face_cells[sub_f_id][1]);
#endif
                }
                else {
#if _DEBUG_
                  bft_printf("bfacecell = %d, _b_face_cells = %d\n",
                             mesh->b_face_cells[f_id], _b_face_cells[sub_f_id]);
#endif
                }
              }
              else {
#if _DEBUG_
                bft_printf("--> sub_face is NOT in cycle, keep standard, "
                           "vtx sign is 0 or +1 :\n");
#endif
                if (is_internal) {
#if _DEBUG_
                  bft_printf("f_id = %d, sub_f_id = %d, "
                             "i_face_cells = %d %d, _i_face_cells = %d %d\n",
                             f_id, sub_f_id,
                             i_face_cells[f_id][0], i_face_cells[f_id][1],
                             _i_face_cells[sub_f_id][0],
                             _i_face_cells[sub_f_id][1]);
#endif
                }
                else {
#if _DEBUG_
                  bft_printf("bfacecell = %d, _b_face_cells = %d\n",
                             mesh->b_face_cells[f_id], _b_face_cells[sub_f_id]);
#endif
                }
              }
            }

          }

        } /* End loop on sub faces */
      } /* End c2f loop */

      /* Convention c_id < c_id_new */
      _i_face_cells[n_i_faces_tot][0] = c_id;
      _i_face_cells[n_i_faces_tot][1] = c_id_new;

      mesh->i_face_vtx_idx[n_i_faces_tot + 1] = n_i_face_vtx_tot + cycle_length;
      n_i_faces_tot++;
      n_i_face_vtx_tot += cycle_length;

      mesh->cell_family[c_id_new] = mesh->cell_family[c_id];

      n_new_cells++;

      n2o_cells[c_id_new] = c_id;

    } /* End loop on cycle */

    /* Free arrays */
    CS_FREE(edge_in_cell);
    CS_FREE(reorder_edge_id);
    CS_FREE(used);
    CS_FREE(cycle_idx);
    CS_FREE(cycle_len);
  } /* End loop on n_cells_cut */
}

/*----------------------------------------------------------------------------*/
/*
 *
 * \brief Rebuild the global interior face connectivity after face subdivision
 *
 * This function replaces the original interior face -> vertex connectivity
 * (i_face_vtx_idx and i_face_vtx_lst) with the new subdivided one.
 *
 * It also rebuilds the face -> cells connectivity (_i_face_cells) by copying
 * the original cell pair for each sub-face generated from the same original face
 *
 * \param[in]      mesh                 pointer to the mesh structure
 * \param[in]      n_i_face_vtx_tot     total number of vertex entries in all
 *                                      sub-faces
 *                                      (size of the final i_face_vtx_lst array)
 * \param[in]      i_face_o2n_idx       old-to-new sub-face index pointer:
 *                                      i_face_o2n_idx[f_id+1] = total number
 *                                      of sub-faces created up to (and including)
 *                                      original face f_id
 * \param[in]      i_sub_face_vtx_idx   sub-face vertex index
 * \param[in]      i_sub_face_vtx_lst   list of vertex ids for sub-faces
 * \param[out]     _i_face_cells        rebuilt face -> cells connectivity array
 *                                      (size = total number of sub-faces)
 */
/*----------------------------------------------------------------------------*/

static void
_update_i_face_connectivity(cs_mesh_t        *mesh,
                            const cs_lnum_t   n_i_face_vtx_tot,
                            const cs_lnum_t   i_face_o2n_idx[],
                            const cs_lnum_t   i_sub_face_vtx_idx[],
                            const cs_lnum_t   i_sub_face_vtx_lst[],
                            cs_lnum_2_t       _i_face_cells[])
{
  CS_FREE(mesh->i_face_vtx_idx);
  CS_FREE(mesh->i_face_vtx_lst);

  const cs_lnum_t n_i_faces_tot = i_face_o2n_idx[mesh->n_i_faces];
  cs_lnum_t *i_face_vtx_idx_glob;
  CS_MALLOC(i_face_vtx_idx_glob, n_i_faces_tot + 1, cs_lnum_t);
  i_face_vtx_idx_glob[0] = 0;

  cs_lnum_t *_i_face_vtx_lst;
  CS_MALLOC(_i_face_vtx_lst, n_i_face_vtx_tot, cs_lnum_t);

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  cs_lnum_t *_i_face_family;
  CS_MALLOC(_i_face_family, n_i_faces_tot, cs_lnum_t);
  const cs_lnum_t *i_face_family = mesh->i_face_family;

  /* Counts to indices */
  cs_lnum_t i_count = 0, i_count_lst = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t s_id = i_face_o2n_idx[f_id];
    cs_lnum_t e_id = i_face_o2n_idx[f_id+1];

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    /* Détail par sous-face locale */
    for (cs_lnum_t s = s_id; s < e_id; s++) {
      cs_lnum_t v0 = i_sub_face_vtx_idx[s];
      cs_lnum_t v1 = i_sub_face_vtx_idx[s+1];
      cs_lnum_t nvtx = v1-v0;

      i_face_vtx_idx_glob[i_count+1] = i_face_vtx_idx_glob[i_count] + nvtx;
      _i_face_cells[i_count][0] = c_id1;
      _i_face_cells[i_count][1] = c_id2;
      _i_face_family[i_count] = i_face_family[f_id];

      i_count++;

      for (cs_lnum_t k = v0; k < v1; k++) {
        _i_face_vtx_lst[i_count_lst++] = i_sub_face_vtx_lst[k];
      }
    }
  }

  assert(i_count == n_i_faces_tot);
  assert(i_count_lst == i_face_vtx_idx_glob[n_i_faces_tot]);
  assert(i_count_lst == n_i_face_vtx_tot);
  mesh->i_face_vtx_idx = i_face_vtx_idx_glob;
  mesh->i_face_vtx_connect_size = i_face_vtx_idx_glob[n_i_faces_tot];
  mesh->i_face_vtx_lst = _i_face_vtx_lst;

  CS_FREE(mesh->i_face_family);
  mesh->i_face_family = _i_face_family;

  /* Update global numbering */

  mesh->n_g_i_faces
    = cs::mesh::o2n_idx_update_global_num(mesh->n_i_faces,
                                          mesh->n_g_i_faces,
                                          i_face_o2n_idx,
                                          &(mesh->global_i_face_num));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Rebuild the global boundary face connectivity after face subdivision
 *
 * This function replaces the original boundary face -> vertex connectivity
 * (b_face_vtx_idx and b_face_vtx_lst) with the new subdivided one
 *
 * It also rebuilds the face -> cell connectivity (_b_face_cells) by copying
 * the original adjacent cell for each sub-face generated from the same
 * original boundary face
 *
 * \param[in]      mesh                 pointer to the mesh structure
 * \param[in]      n_vtx_tot            total number of vertex entries in all
 *                                      sub-faces
 *                                      (size of the final b_face_vtx_lst array)
 * \param[in]      b_face_o2n_idx       old-to-new sub-face index pointer:
 *                                      b_face_o2n_idx[f_id+1] = total number
 *                                      of sub-faces created up to (and including)
 *                                      original face f_id
 * \param[in]      b_sub_face_vtx_idx   vertex index of sub-face
 * \param[in]      b_sub_face_vtx_lst   vertex ids for sub-faces
 * \param[out]     _b_face_cells        rebuilt face -> cell connectivity array
 *                                      (size = total number of sub-faces)
 */
/*----------------------------------------------------------------------------*/

static void
_update_b_face_connectivity(cs_mesh_t        *mesh,
                            const cs_lnum_t   n_vtx_tot,
                            const cs_lnum_t   b_face_o2n_idx[],
                            const cs_lnum_t   b_sub_face_vtx_idx[],
                            const cs_lnum_t   b_sub_face_vtx_lst[],
                            cs_lnum_t         _b_face_cells[])
{
  CS_FREE(mesh->b_face_vtx_idx);
  CS_FREE(mesh->b_face_vtx_lst);

  const cs_lnum_t n_b_faces_tot = b_face_o2n_idx[mesh->n_b_faces];
  cs_lnum_t *b_face_vtx_idx_glob;
  CS_MALLOC(b_face_vtx_idx_glob, n_b_faces_tot + 1, cs_lnum_t);
  b_face_vtx_idx_glob[0] = 0;

  cs_lnum_t *_b_face_vtx_lst = nullptr;
  CS_MALLOC(_b_face_vtx_lst, n_vtx_tot, cs_lnum_t);

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  cs_lnum_t *_b_face_family;
  CS_MALLOC(_b_face_family, n_b_faces_tot, cs_lnum_t);
  const cs_lnum_t *b_face_family = mesh->b_face_family;

  /* Counts to indices */
  cs_lnum_t count = 0, count_lst = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t s_id = b_face_o2n_idx[f_id];
    cs_lnum_t e_id = b_face_o2n_idx[f_id+1];

    cs_lnum_t c_id = b_face_cells[f_id];

    /* Local sub-faces details */
    for (cs_lnum_t s = s_id; s < e_id; s++) {
      cs_lnum_t v0 = b_sub_face_vtx_idx[s];
      cs_lnum_t v1 = b_sub_face_vtx_idx[s+1];
      cs_lnum_t nvtx = v1-v0;

      b_face_vtx_idx_glob[count+1] = b_face_vtx_idx_glob[count] + nvtx;
      _b_face_cells[count] = c_id;
      _b_face_family[count] = b_face_family[f_id];

      count++;

      for (cs_lnum_t k = v0; k < v1; k++) {
        _b_face_vtx_lst[count_lst++] = b_sub_face_vtx_lst[k];
      }
    }
  }

  assert(count == n_b_faces_tot);
  assert(count_lst == b_face_vtx_idx_glob[n_b_faces_tot]);
  assert(count_lst == n_vtx_tot);
  mesh->b_face_vtx_idx = b_face_vtx_idx_glob;
  mesh->b_face_vtx_connect_size = b_face_vtx_idx_glob[n_b_faces_tot];
  mesh->b_face_vtx_lst = _b_face_vtx_lst;

  CS_REALLOC(mesh->b_face_family, n_b_faces_tot, cs_lnum_t);
  CS_REALLOC(mesh->b_face_r_c_idx, n_b_faces_tot, char);
  for (cs_lnum_t i = mesh->n_b_faces; i < n_b_faces_tot; i++) {
    mesh->b_face_family[i] = 1;
    if (mesh->have_r_gen)
      mesh->b_face_r_c_idx[i] = 127;
  }

  CS_FREE(mesh->b_face_family);
  mesh->b_face_family = _b_face_family;

  mesh->n_g_b_faces
    = cs::mesh::o2n_idx_update_global_num(mesh->n_b_faces,
                                          mesh->n_g_b_faces,
                                          b_face_o2n_idx,
                                          &(mesh->global_b_face_num));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Cut (subdivide) a single face intersected by the cutting plane(s).
 *
 * This function detects if the face is cut (i.e. has valid intersection
 * points on its edges), and if so, subdivides it into one or more sub-faces
 * (polygons) that lie entirely on one side or the other of the cutting plane.
 *
 * It also identifies and records "dark" and "light" edges:
 *   - dark  edges: edges connecting a negative side vertex to a cut edge
 *   - light edges: edges connecting a positive side vertex to a cut edge
 *
 * These dark/light edges are used later to help close polygons and detect
 * intersection contours (typically for immersed boundary or cut-cell logic).
 *
 * Special cases:
 *   - No intersection -> the original face is kept unchanged (1 sub-face)
 *   - Intersections present -> the face is split into multiple sub-faces
 *
 * \param[in]      f_id                     face id
 * \param[in]      f2e_idx                  face to edge index (start pointers)
 * \param[in]      f2e_ids                  face to edge connectivity list
 * \param[in]      f2v_idx                  face to vertex index
 * \param[in]      f2v_ids                  face to vertex connectivity list
 * \param[in]      v2v_ids                  vertex to vertex connectivity list
 * \param[in]      edge_v0                  for each edge, index of vertex with
 *                                          smaller id
 * \param[in]      e_v_idx                  start index of new vertices per edge
 *                                          (e_v_idx[edge_id+1] > e_v_idx[edge_id]
 *                                            --> cut)
 * \param[in]      vertex_sign              sign of each vertex
 * \param[in]      intx                     intersection parameter t \in [0,1]
 *                                          on each edge (-1 if not intersected)
 * \param[out]     face_o2n_idx             new (sub-face) index pointer for
 *                                          this face
 *                                          (cumulative: face_o2n_idx[f_id+1] =
 *                                          total sub-faces so far)
 * \param[out]     face_o2n_connect_idx     cumulative vertex connectivity
 *                                          index after this face
 * \param[in,out]  n_vtx_tot                total number of vertices written in
 *                                          sub_face_vtx_lst so far
 * \param[in,out]  n_sub_face_tot           total number of sub-faces created
 * \param[out]     sub_face_vtx_idx         so far sub-face to vertex start
 *                                          indices (compact format)
 * \param[out]     sub_face_vtx_lst         concatenated list of all
 *                                          sub-face vertex ids
 * \param[in,out]  n_dark_edge_in_face_tot  total number of dark edges found
 * \param[out]     dark_edge_in_face_idx    so far cumulative index of dark
 *                                          edges per original face
 * \param[out]     dark_edge_in_face        array of [edge_id_a, edge_id_b] pairs
 *                                          for dark edges
 * \param[in,out]  n_light_edge_in_face_tot total number of light edges found
 * \param[out]     light_edge_in_face_idx   so far cumulative index of light
 *                                          adges per original face
 * \param[out]     light_edge_in_face       array of [edge_id_a, edge_id_b] pairs
 *                                          for light edges
 */
/*----------------------------------------------------------------------------*/

static void
_cut_face(const cs_lnum_t  f_id,
          const cs_lnum_t  f2e_idx[],
          const cs_lnum_t  f2e_ids[],
          const cs_lnum_t  f2v_idx[],
          const cs_lnum_t  f2v_ids[],
          const cs_lnum_t  v2v_ids[],
          const cs_lnum_t  edge_v0[],
          const cs_lnum_t  e_v_idx[],
          const int        vertex_sign[],
          const double     intx[],
          cs_lnum_t        face_o2n_idx[],
          cs_lnum_t        face_o2n_connect_idx[],
          cs_lnum_t       &n_vtx_tot,
          cs_lnum_t       &n_sub_face_tot,
          cs_lnum_t        sub_face_vtx_idx[],
          cs_lnum_t        sub_face_vtx_lst[],
          cs_lnum_t       &n_dark_edge_in_face_tot,
          cs_lnum_t        dark_edge_in_face_idx[],
          cs_lnum_t        dark_edge_in_face[][2],
          cs_lnum_t       &n_light_edge_in_face_tot,
          cs_lnum_t        light_edge_in_face_idx[],
          cs_lnum_t        light_edge_in_face[][2])
{
  int n_intx_in_face = 0;
  int n_dark_edge_in_face = 0, n_light_edge_in_face = 0;
  int n_sub_face = 0;

  const cs_lnum_t s_id_f = f2e_idx[f_id];
  const cs_lnum_t e_id_f = f2e_idx[f_id + 1];

  for (cs_lnum_t fidx = s_id_f; fidx < e_id_f; fidx++) {

    cs_lnum_t edge_id = f2e_ids[fidx];
    cs_real_t _intx = intx[edge_id];
    if (_intx >= 0 && _intx <= 1)
      n_intx_in_face++;
  }

  if (n_intx_in_face > 0)
    assert(n_intx_in_face > 1);

  //const cs_lnum_t max_sub_edge_in_face = 4*n_intx_in_face;
  cs_lnum_t edge_intx_vfluid[100][2]; //FIXME
  int n_fluid_edge_in_face = 0;

  if (n_intx_in_face == 0) {
    n_sub_face++;
    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t n_vtx = f2v_idx[f_id+1] - s_id;

    for (cs_lnum_t i = 0; i < n_vtx; ++i)
      sub_face_vtx_lst[n_vtx_tot+i] = f2v_ids[s_id + i];

    n_sub_face_tot += n_sub_face;
    n_vtx_tot += n_vtx;

    sub_face_vtx_idx[n_sub_face_tot] = n_vtx_tot;
    face_o2n_idx[f_id+1] = n_sub_face_tot;
    face_o2n_connect_idx[f_id+1] = n_vtx_tot;
    dark_edge_in_face_idx[f_id+1] = n_dark_edge_in_face_tot;
    light_edge_in_face_idx[f_id+1] = n_light_edge_in_face_tot;
    return; // Go to the next face_id
  }

  // here n_intx_in_face > 2 : Need to create sub_face and vertex contour

  for (cs_lnum_t fidx = s_id_f; fidx < e_id_f; fidx++) {

    cs_lnum_t edge_id = f2e_ids[fidx];
    cs_real_t _intx = intx[edge_id];

    if (_intx >= 0 && _intx <= 1) {

      cs_lnum_t n_vtx = 0;

      sub_face_vtx_lst[n_vtx_tot + n_vtx] = e_v_idx[edge_id];
      n_vtx++;

      sub_face_vtx_idx[n_sub_face_tot + n_sub_face] = n_vtx_tot;

      cs_lnum_t v0 = edge_v0[edge_id];
      cs_lnum_t v1 = v2v_ids[edge_id];

      cs_lnum_t _fidx = fidx+1;
      if (_fidx >= e_id_f)
        _fidx -= (e_id_f - s_id_f);

      cs_lnum_t fidx_start = _fidx;

      bool init_fluid_edge = false;

      /* Loop over following edges */
      for (cs_lnum_t fidx_1 = fidx_start;
          fidx_1 < fidx_start + (e_id_f - s_id_f) - 1; fidx_1++) {

        cs_lnum_t _fidx_1 = fidx_1;
        if (_fidx_1 >= e_id_f)
          _fidx_1 -= (e_id_f - s_id_f);

        cs_lnum_t edge_id_1 = f2e_ids[_fidx_1];
        cs_lnum_t v0_1 = edge_v0[edge_id_1];
        cs_lnum_t v1_1 = v2v_ids[edge_id_1];

        cs_lnum_t common_vi = -1;
        if (v0 == v0_1 || v0 == v1_1)
          common_vi = v0;
        if (v1 == v0_1 || v1 == v1_1)
          common_vi = v1;

        sub_face_vtx_lst[n_vtx_tot + n_vtx] = common_vi;
        n_vtx++;

        /* Process dark edge */
        assert(common_vi > -1);
        if (vertex_sign[common_vi] == VTX_FLUID) {

          /* If the first common_vi is fluid -> not construct sub faces */
          for (cs_lnum_t i = 0; i < n_vtx; i++) {
            //sub_face_vtx_lst[n_vtx_tot + n_vtx] = -1;
            sub_face_vtx_lst[n_vtx_tot + i] = -1;
          }
          n_vtx = 0;

          /* Save fluid internal edges connected to fluid common vertex */
          if (!init_fluid_edge) {
            edge_intx_vfluid[n_fluid_edge_in_face][0] = e_v_idx[edge_id];
            edge_intx_vfluid[n_fluid_edge_in_face][1] = common_vi;
            n_fluid_edge_in_face++;
            init_fluid_edge = true;
          }

          cs_real_t _intx_1 = intx[edge_id_1];

          if (!(_intx_1 >= 0 && _intx_1 <= 1)) {

            /* Save fluid internal edges connected to fluid common vertex */
            edge_intx_vfluid[n_fluid_edge_in_face][0] = v0_1;
            edge_intx_vfluid[n_fluid_edge_in_face][1] = v1_1;

            v0 = v0_1;
            v1 = v1_1;
            n_fluid_edge_in_face++;
            assert(n_fluid_edge_in_face <= 100);

            continue;
          }
          else {

            /* Save fluid internal edges connected to fluid common vertex */
            edge_intx_vfluid[n_fluid_edge_in_face][0] = common_vi;
            edge_intx_vfluid[n_fluid_edge_in_face][1] = e_v_idx[edge_id_1];
            n_fluid_edge_in_face++;

            cs_lnum_t id = n_dark_edge_in_face_tot + n_dark_edge_in_face;
            dark_edge_in_face[id][0] = edge_id_1;
            dark_edge_in_face[id][1] = edge_id;

            n_dark_edge_in_face++;
            assert(n_fluid_edge_in_face <= 100);

            break;
          }
        }

        /* Process light edge (vertex_sign[common_vi] = 1) */

        assert(vertex_sign[common_vi] == VTX_SOLID);

        cs_real_t _intx_1 = intx[edge_id_1];

        if (!(_intx_1 >= 0 && _intx_1 <= 1)) {
          v0 = v0_1;
          v1 = v1_1;
          continue;
        }
        else {
          /* edge creation */
          cs_lnum_t id = n_light_edge_in_face_tot + n_light_edge_in_face;

          light_edge_in_face[id][0] = edge_id;
          light_edge_in_face[id][1] = edge_id_1;

          n_light_edge_in_face++;

          sub_face_vtx_lst[n_vtx_tot + n_vtx] = e_v_idx[edge_id_1];
          n_vtx++;
          n_vtx_tot += n_vtx;

          sub_face_vtx_idx[n_sub_face_tot + n_sub_face + 1] = n_vtx_tot;
          n_sub_face++;

          break;
        }
      }

    } /* End test on intx */
    else {
#if _DEBUG_
      bft_printf("face_id = %d, edge_id = %d : NO_INTX --> continue\n",
                 f_id, edge_id);
#endif
    }

  } /* End loop face to edge */

  n_dark_edge_in_face_tot += n_dark_edge_in_face;
  n_light_edge_in_face_tot += n_light_edge_in_face;
  dark_edge_in_face_idx[f_id+1] = n_dark_edge_in_face_tot;
  light_edge_in_face_idx[f_id+1] = n_light_edge_in_face_tot;

  /* Contruct fluid contour by connecting light edges with fluid edges */

  const cs_lnum_t le0 = light_edge_in_face_idx[f_id];
  const cs_lnum_t le1 = light_edge_in_face_idx[f_id+1];
  const cs_lnum_t n_le = le1 - le0;

  // Add light edge to fluid edges and reordering
  for (cs_lnum_t i = 0; i < n_le; i++) {
    edge_intx_vfluid[i + n_fluid_edge_in_face][0]
      = e_v_idx[light_edge_in_face[i+le0][0]];
    edge_intx_vfluid[i + n_fluid_edge_in_face][1]
      = e_v_idx[light_edge_in_face[i+le0][1]];
  }
  n_fluid_edge_in_face += n_light_edge_in_face;

  cs_lnum_t count = 0;
  // Take first fluid edge as reference
  for (count = 0; count < 2; count++)
    sub_face_vtx_lst[n_vtx_tot + count] = edge_intx_vfluid[0][count];

  while (count != n_fluid_edge_in_face) {
    for (cs_lnum_t j = 1; j < n_fluid_edge_in_face; j++) {

      for (cs_lnum_t k = 0; k < 2; k++) {
        if (edge_intx_vfluid[j][k] == sub_face_vtx_lst[n_vtx_tot + count - 1]) {
          sub_face_vtx_lst[n_vtx_tot + count]
            = edge_intx_vfluid[j][(k+1)%2];
          count++;
          break;
        }
      }

      if (count == n_fluid_edge_in_face)
        break;
    }
  }

  sub_face_vtx_idx[n_sub_face_tot + n_sub_face + 1] = n_vtx_tot + count;
  n_sub_face++;
  n_vtx_tot += count;

  n_sub_face_tot += n_sub_face;
  face_o2n_idx[f_id+1] = n_sub_face_tot;
  face_o2n_connect_idx[f_id+1] = n_vtx_tot;

  /* Debug print */

#if _DEBUG_

  // Global sub-faces for face f_id : [s_id, e_id)
  cs_lnum_t s_id = face_o2n_idx[f_id];
  cs_lnum_t e_id = face_o2n_idx[f_id+1];

  // Global range connectivity for face f_id
  cs_lnum_t v_beg = sub_face_vtx_idx[s_id];
  cs_lnum_t v_end = sub_face_vtx_idx[e_id];

  cs_lnum_t n_sf = e_id - s_id;
  cs_lnum_t n_vtx_face = v_end - v_beg;

  printf("[face %d] subfaces: s_id=%d, e_id=%d (count=%d), "
         "connect range: v_id=%d..%d (count=%d)\n",
         f_id, s_id, e_id, n_sf,
         v_beg, (v_end > 0 ? v_end - 1 : 0), n_vtx_face);

  // Local sub face details
  for (cs_lnum_t s = s_id; s < e_id; ++s) {
    cs_lnum_t v0 = sub_face_vtx_idx[s];
    cs_lnum_t v1 = sub_face_vtx_idx[s+1];
    cs_lnum_t sf_local = s - s_id;

    printf("  subface_local=%d (global=%d): vtx s_id=%d, e_id=%d -> ",
           sf_local, s, v0, v1);

    for (cs_lnum_t k = v0; k < v1; ++k)
      printf("%d ", sub_face_vtx_lst[k]);
    printf("\n");
  }

  // Dark edges for the face

  const cs_lnum_t de0 = dark_edge_in_face_idx[f_id];
  const cs_lnum_t de1 = dark_edge_in_face_idx[f_id+1];
  const cs_lnum_t n_de = de1 - de0;

  if (de1 > de0) {
    printf("  dark_edges: s_id=%d, e_id=%d (count=%d): ",
           de0, de1, n_de);
    for (cs_lnum_t i = de0; i < de1; ++i)
      printf("(%d,%d) ", dark_edge_in_face[i][0], dark_edge_in_face[i][1]);
    printf("\n");
  }

  // Light edges for the face

  if (le1 > le0) {
    printf("  light_edges: s_id=%d, e_id=%d (count=%d): ",
           le0, le1, n_le);
    for (cs_lnum_t i = le0; i < le1; ++i)
      printf("(%d,%d) ", light_edge_in_face[i][0], light_edge_in_face[i][1]);
    printf("\n");
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Subdivide all faces (boundary and interior) of cells intersected by
 *        the cutting plane(s).
 *
 *    - Calls _cut_face() on every boundary face and every interior
 *      face of the mesh.
 *    - For each face:
 *      - If no valid intersection -> keeps the original face (1 sub-face)
 *      - If intersections exist -> splits the face into multiple sub-faces
 *      - Builds connectivity lists for sub-faces (vertex indices and lists)
 *      - Collects and stores "light" and "dark" edges per original face
 *        (used later to close intersection polygons and assign sub-faces
 *         to new cells)
 *    - Accumulates total number of sub-faces and vertex entries for boundary
 *      and interior faces.
 *
 * The output arrays (o2n_idx, sub_face_vtx_idx/lst, light/dark edge arrays) are
 * filled in compact/indexed format and serve as input for later steps:
 *   - connectivity rebuild
 * (_update_b_face_connectivity / _update_i_face_connectivity)
 *   - new cell and immersed face creation (_cut_cells)
 *
 * \param[in,out] mesh                      pointer to the mesh structure
 * \param[in]     intx                      intersection parameter t ∈ [0,1]
 *                                          on each edge (-1 if not cut)
 * \param[in]     v2v                       vertex to vertex adjacency
 * \param[in]     f2e                       face to edge adjacency
 * \param[in]     edge_v0                   for each edge, index of vertex with
 *                                          smaller id
 * \param[in,out] vertex_sign               vertex sign classification
 * \param[in]     cell_flag                 flag array marking cells disabled
 *                                          due to inconsistent vertex signs
 * \param[in]     e_v_idx                   global start id of new vertices
 *                                          created on each edge
 * \param[in,out] n_sub_b_face_tot          total number of boundary sub-faces
 *                                          created
 * \param[in,out] n_b_face_vtx_tot          total number of vertex entries in
 *                                          boundary sub-faces
 * \param[in,out] n_sub_i_face_tot          total number of interior sub-faces
*                                           created
 * \param[in,out] n_i_face_vtx_tot          total number of vertex entries in
 *                                          interior sub-faces
 * \param[out]    b_face_o2n_idx            cumulative sub-face count per
 *                                          original boundary face
 * \param[out]    b_face_o2n_connect_idx    cumulative vertex count per original
 *                                          boundary face
 * \param[out]    i_face_o2n_idx            cumulative sub-face count per
 *                                          original interior face
 * \param[out]    i_face_o2n_connect_idx    cumulative vertex count per original
 *                                          interior face
 * \param[out]    b_sub_face_vtx_idx        vertex start indices for boundary
 *                                          sub-faces
 * \param[out]    b_sub_face_vtx_lst        concatenated vertex list for boundary
 *                                          sub-faces
 * \param[out]    i_sub_face_vtx_idx        vertex start indices for interior
 *                                          sub-faces
 * \param[out]    i_sub_face_vtx_lst        concatenated vertex list for interior
 *                                          sub-faces
 * \param[out]    light_edge_in_b_face_idx  cumulative start index of light edges
 *                                          per boundary face
 * \param[out]    light_edge_in_i_face_idx  cumulative start index of light edges
 *                                          per interior face
 * \param[out]    light_edge_in_b_face      [edge_a, edge_b] pairs for light
 *                                           edges on boundary faces
 * \param[out]    light_edge_in_i_face      [edge_a, edge_b] pairs for light
 *                                           edges on interior faces
 * \param[out]    dark_edge_in_b_face_idx   cumulative start index of dark edges
 *                                          per boundary face
 * \param[out]    dark_edge_in_i_face_idx   cumulative start index of dark edges
 *                                          per interior face
 * \param[out]    dark_edge_in_b_face       [edge_a, edge_b] pairs for dark edges
 *                                           on boundary faces
 * \param[out]    dark_edge_in_i_face       [edge_a, edge_b] pairs for dark edges
 *                                           on interior faces
 */
/*----------------------------------------------------------------------------*/

static void
_cut_faces(cs_mesh_t            *mesh,
           const cs_real_t       intx[],
           const cs_adjacency_t *v2v,
           const cs_adjacency_t *f2e,
           const cs_lnum_t       edge_v0[],
           const int             vertex_sign[],
           const int             cell_flag[],
           const cs_lnum_t       e_v_idx[],
           cs_lnum_t            &n_sub_b_face_tot,
           cs_lnum_t            &n_b_face_vtx_tot,
           cs_lnum_t            &n_sub_i_face_tot,
           cs_lnum_t            &n_i_face_vtx_tot,
           cs_lnum_t             b_face_o2n_idx[],
           cs_lnum_t             b_face_o2n_connect_idx[],
           cs_lnum_t             i_face_o2n_idx[],
           cs_lnum_t             i_face_o2n_connect_idx[],
           cs_lnum_t             b_sub_face_vtx_idx[],
           cs_lnum_t             b_sub_face_vtx_lst[],
           cs_lnum_t             i_sub_face_vtx_idx[],
           cs_lnum_t             i_sub_face_vtx_lst[],
           cs_lnum_t             light_edge_in_b_face_idx[],
           cs_lnum_t             light_edge_in_i_face_idx[],
           cs_lnum_2_t           light_edge_in_b_face[],
           cs_lnum_2_t           light_edge_in_i_face[],
           cs_lnum_t             dark_edge_in_b_face_idx[],
           cs_lnum_t             dark_edge_in_i_face_idx[],
           cs_lnum_2_t           dark_edge_in_b_face[],
           cs_lnum_2_t           dark_edge_in_i_face[])
{
  cs_lnum_t *v2v_ids = v2v->ids;

  cs_lnum_t n_dark_edge_in_b_face_tot = 0;
  cs_lnum_t n_light_edge_in_b_face_tot = 0;

  cs_lnum_t *b_f2v_idx = mesh->b_face_vtx_idx;
  cs_lnum_t *b_f2v_ids = mesh->b_face_vtx_lst;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {

    cs_lnum_t c_id = mesh->b_face_cells[f_id];
    if (cell_flag[c_id] == 1) {

      /* Do not cut faces shared by disabled cells */

      cs_lnum_t n_sub_face = 1;
      cs_lnum_t s_id = b_f2v_idx[f_id];
      cs_lnum_t n_vtx = b_f2v_idx[f_id+1] - s_id;

      for (cs_lnum_t i = 0; i < n_vtx; ++i)
        b_sub_face_vtx_lst[n_b_face_vtx_tot+i] = b_f2v_ids[s_id + i];

      n_sub_b_face_tot += n_sub_face;
      n_b_face_vtx_tot += n_vtx;

      b_sub_face_vtx_idx[n_sub_b_face_tot] = n_b_face_vtx_tot;
      b_face_o2n_idx[f_id+1] = n_sub_b_face_tot;
      b_face_o2n_connect_idx[f_id+1] = n_b_face_vtx_tot;
      dark_edge_in_b_face_idx[f_id+1] = n_dark_edge_in_b_face_tot;
      light_edge_in_b_face_idx[f_id+1] = n_light_edge_in_b_face_tot;
      continue;
    }

    _cut_face(f_id,
              f2e->idx + mesh->n_i_faces,
              f2e->ids,
              b_f2v_idx,
              b_f2v_ids,
              v2v_ids,
              edge_v0,
              e_v_idx,
              vertex_sign,
              intx,
              b_face_o2n_idx,
              b_face_o2n_connect_idx,
              n_b_face_vtx_tot,
              n_sub_b_face_tot,
              b_sub_face_vtx_idx,
              b_sub_face_vtx_lst,
              n_dark_edge_in_b_face_tot,
              dark_edge_in_b_face_idx,
              dark_edge_in_b_face,
              n_light_edge_in_b_face_tot,
              light_edge_in_b_face_idx,
              light_edge_in_b_face);

  } /* End boundary face loop */

  cs_lnum_t n_dark_edge_in_i_face_tot = 0;
  cs_lnum_t n_light_edge_in_i_face_tot = 0;

  cs_lnum_t *i_f2v_idx = mesh->i_face_vtx_idx;
  cs_lnum_t *i_f2v_ids = mesh->i_face_vtx_lst;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {

    cs_lnum_t c_id0 = mesh->i_face_cells[f_id][0];
    cs_lnum_t c_id1 = mesh->i_face_cells[f_id][1];

    if (cell_flag[c_id0] == 1 && cell_flag[c_id1] == 1) {

      /* Do not cut faces shared by disabled cells */

      cs_lnum_t n_sub_face = 1;
      cs_lnum_t s_id = i_f2v_idx[f_id];
      cs_lnum_t n_vtx = i_f2v_idx[f_id+1] - s_id;

      for (cs_lnum_t i = 0; i < n_vtx; ++i)
        i_sub_face_vtx_lst[n_i_face_vtx_tot+i] = i_f2v_ids[s_id + i];

      n_sub_i_face_tot += n_sub_face;
      n_i_face_vtx_tot += n_vtx;

      i_sub_face_vtx_idx[n_sub_i_face_tot] = n_i_face_vtx_tot;
      i_face_o2n_idx[f_id+1] = n_sub_i_face_tot;
      i_face_o2n_connect_idx[f_id+1] = n_i_face_vtx_tot;
      dark_edge_in_i_face_idx[f_id+1] = n_dark_edge_in_i_face_tot;
      light_edge_in_i_face_idx[f_id+1] = n_light_edge_in_i_face_tot;
      continue;
    }

    _cut_face(f_id,
              f2e->idx,
              f2e->ids,
              i_f2v_idx,
              i_f2v_ids,
              v2v_ids,
              edge_v0,
              e_v_idx,
              vertex_sign,
              intx,
              i_face_o2n_idx,
              i_face_o2n_connect_idx,
              n_i_face_vtx_tot,
              n_sub_i_face_tot,
              i_sub_face_vtx_idx,
              i_sub_face_vtx_lst,
              n_dark_edge_in_i_face_tot,
              dark_edge_in_i_face_idx,
              dark_edge_in_i_face,
              n_light_edge_in_i_face_tot,
              light_edge_in_i_face_idx,
              light_edge_in_i_face);

  } /* End internal face loop */

  // init i_face_cells in the halo for cs_mesh_init_halo
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (mesh->i_face_cells[f_id][j] >= mesh->n_cells) {
        mesh->i_face_cells[f_id][j] = -1;
        break;
      }
    }
  }

  bft_printf(" Boundary and internal faces successfully subdivided\n");

} /* End loop on n_cells_cut */

/*----------------------------------------------------------------------------*/
/*
 * \brief Build new vertices on selected edges.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The coordinates and numbering arrays should be resized before calling
 * this function (to allow for vertices inserted on edges, faces, and
 * cells with a single resize).
 *
 * \param[in]  m          mesh
 * \param[in]  v2v        vertex->vertex adjacency
 * \param[in]  n_e_vtx    local number of vertices added on edges
 * \param[in]  e_v_idx    for each edge, start index of added vertices
 * \param[in]  g_e_v_num  for each edge, start global number of
 *                        added vertices (1 to n, to shift) or nullptr
 *
 * \return  local number of vertices added on mid edges
 */
/*----------------------------------------------------------------------------*/

static void
_build_edge_vertices(cs_mesh_t             *m,
                     const cs_adjacency_t  *v2v,
                     double                 intx[],
                     const cs_lnum_t        n_e_vtx,
                     const cs_lnum_t        e_v_idx[],
                     const cs_gnum_t       *g_e_v_num)
{
  const cs_lnum_t n_edges = v2v->idx[v2v->n_elts];
  cs_real_3_t *vtx_coord = (cs_real_3_t *)m->vtx_coord;

  if (m->global_vtx_num != nullptr && g_e_v_num != nullptr) {

    cs_gnum_t gnum_shift = m->n_g_vertices - 1;

    for (cs_lnum_t id0 = 0; id0 < v2v->n_elts; id0++) {
      for (cs_lnum_t j = v2v->idx[id0]; j < v2v->idx[id0+1]; j++) {
        cs_lnum_t s_id = e_v_idx[j];
        cs_lnum_t e_id = e_v_idx[j+1];
        cs_lnum_t n_sub = e_id - s_id;
        assert(n_sub == 0 || n_sub == 1);
        double _intx = intx[j];
        cs_lnum_t id1 = v2v->ids[j];

        for (cs_lnum_t k = 0; k < n_sub; k++) {
          cs_lnum_t id2 = k + s_id;

          for (cs_lnum_t l = 0; l < 3; l++)
            vtx_coord[id2][l]
              =  vtx_coord[id0][l]
               + _intx * (vtx_coord[id1][l] - vtx_coord[id0][l]);
          m->global_vtx_num[id2+k] = gnum_shift + g_e_v_num[j] + (cs_gnum_t)k;
        }
      }
    }

  }
  else {

    for (cs_lnum_t id0 = 0; id0 < v2v->n_elts; id0++) {
      for (cs_lnum_t j = v2v->idx[id0]; j < v2v->idx[id0+1]; j++) {
        cs_lnum_t s_id = e_v_idx[j];
        cs_lnum_t e_id = e_v_idx[j+1];
        cs_lnum_t n_sub = e_id - s_id;
        assert(n_sub == 0 || n_sub == 1);
        double _intx = intx[j];
        cs_lnum_t id1 = v2v->ids[j];
        for (cs_lnum_t k = 0; k < n_sub; k++) {
          cs_lnum_t id2 = k + s_id;
          for (cs_lnum_t l = 0; l < 3; l++)
            vtx_coord[id2][l]
              =  vtx_coord[id0][l]
               + _intx * (vtx_coord[id1][l] - vtx_coord[id0][l]);
        }
      }
    }

    if (m->global_vtx_num != nullptr) {
      for (cs_lnum_t i = 0; i < n_e_vtx; i++)
        m->global_vtx_num[m->n_vertices + i] = m->n_vertices + i + 1;
    }

  }

  m->n_vertices = cs::max(m->n_vertices, e_v_idx[n_edges]);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Function to compute the normal of a triangle defined by 3 points
 *
 * parameters:
 *   coords <-- triangle coordinates
 *   normal --> triangle normal
*/
/*----------------------------------------------------------------------------*/

static void
_compute_normal(cs_real_t  coords[3][3],
                cs_real_t  normal[3])
{
  /* Compute and Write the face normal coordinates */
  const cs_real_t *a = coords[0], *b = coords[1], *c = coords[2];

  // ab vect ac
  normal[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
  normal[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
  normal[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

  cs_real_t norm = cs_math_3_norm(normal);

  for (int dir = 0; dir < 3; dir ++)
    normal[dir] /= norm;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute intersection between edges of the main mesh and STL triangles
 *
 * \param[in]      stl_mesh     pointer to STL mesh
 * \param[in,out]  mesh         pointer to the mesh structure to modify
 * \param[in]      v2v          vertex to vertex connectivity
 * \param[in]      c2e          cell to edge connectivity
 * \param[in]      edge_v0      for each edge_id, index of vertex v0
 * \param[in,out]  intx         intersection parameter t \in [0,1] on each
 *                              edge (-1 if not intersected)
 * \param[in,out]  vertex_sign  vertex sign (fluid, solid, interface)
 * \param[in,out]  e_v_idx      for each edge, start index of added vertices
 * \param[in,out]  n_c_cut      number of cells detected as intersected
 *                              (≥ 3 intersections)
 * \param[in,out]  c_cut        array of cell ids that are cut
 *                              (size ≥ n_c_cut)
 * \param[in,out]  cell_flag    flag array marking cells disabled due to
 *                              inconsistent vertex signs
 *
 * Return: the vertex_sign array reallocated with the added vertex and
 *         propagated only in the cutted cells
 */
/*----------------------------------------------------------------------------*/

static int *
_cut_edges(const cs_stl_mesh_t    *stl_mesh,
           cs_mesh_t              *mesh,
           const cs_adjacency_t   *v2v,
           const cs_adjacency_t   *c2e,
           cs_lnum_t              edge_v0[],
           cs_real_t              intx[],
           int                    vertex_sign[],
           cs_lnum_t              e_v_idx[],
           cs_lnum_t             &n_c_cut,
           cs_lnum_t              c_cut[],
           int                   *cell_flag)
{
  const cs_lnum_t n_vtx = mesh->n_vertices;
  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)mesh->vtx_coord;
  cs_lnum_t *v2v_idx = v2v->idx;
  cs_lnum_t *v2v_ids = v2v->ids;
  cs_lnum_t n_edges = v2v_idx[n_vtx];

  cs_real_t *best_angle;
  CS_MALLOC(best_angle, n_vtx, cs_real_t);

  for (cs_lnum_t v0 = 0; v0 < n_vtx; v0++) {
    cs_lnum_t s_id = v2v_idx[v0];
    cs_lnum_t e_id = v2v_idx[v0+1];
    vertex_sign[v0] = VTX_UNKNOWN;
    best_angle[v0] = -1.0;

    for (cs_lnum_t vidx = s_id; vidx < e_id; vidx++) {
      edge_v0[vidx] = v0;
      intx[vidx] = -1.0;
      e_v_idx[vidx+1] = 0;
    }
  }

  const cs_real_3_t *stl_vtx_coord = stl_mesh->coords;
  const cs_lnum_3_t *tria_vtx_ids = stl_mesh->tria_vtx_ids;

  cs_lnum_t n_input_edges = n_edges;
  cs_lnum_t n_selected_edges = 0;
  cs_lnum_t *input_edges = nullptr;
  cs_lnum_t *selected_edges = nullptr;
  cs_lnum_t *tria_in_edges_lst = nullptr;
  cs_lnum_t *tria_in_edges_idx = nullptr;
  cs_lnum_t *edges_selected_idx = nullptr;

  cs_lnum_t max_size = stl_mesh->n_faces;

  CS_MALLOC(input_edges, n_edges, cs_lnum_t);
  CS_MALLOC(selected_edges, n_edges, cs_lnum_t);
  CS_MALLOC(tria_in_edges_idx, n_edges, cs_lnum_t);
  CS_MALLOC(tria_in_edges_lst, max_size, cs_lnum_t);
  CS_MALLOC(edges_selected_idx, n_edges, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    input_edges[i] = i;
    edges_selected_idx[i] = -1;
  }

  cs_stl_intersection(stl_mesh,
                      CS_MESH_LOCATION_VERTICES,
                      v2v,
                      n_input_edges,
                      input_edges,
                      &n_selected_edges,
                      selected_edges,
                      tria_in_edges_idx,
                      &tria_in_edges_lst,
                      &max_size);

  CS_FREE(input_edges);
  CS_FREE(edges_selected_idx);

  bft_printf(" List of intersected edges completed\n");
  bft_printf("   Total number of edges in the mesh: %d\n", n_edges);
  bft_printf("   Number of intersected edges: %d\n\n", n_selected_edges);

  /* compute normals of STL triangles */
  cs_real_t *stl_normals;
  CS_MALLOC(stl_normals, stl_mesh->n_faces*3, cs_real_t);

  for (cs_lnum_t i = 0; i < stl_mesh->n_faces; i++) {
    cs_real_t coords_loc[3][3];

    const cs_lnum_t vtx_ids_l[3] = {tria_vtx_ids[i][0],
                                    tria_vtx_ids[i][1],
                                    tria_vtx_ids[i][2]};

    for (int dir = 0; dir < 3; dir ++) {
      coords_loc[0][dir] = stl_vtx_coord[vtx_ids_l[0]][dir];
      coords_loc[1][dir] = stl_vtx_coord[vtx_ids_l[1]][dir];
      coords_loc[2][dir] = stl_vtx_coord[vtx_ids_l[2]][dir];
    }

    _compute_normal(coords_loc, stl_normals + (3*i));
  }

  for (cs_lnum_t e_id = 0; e_id < n_selected_edges; e_id++) {
    const cs_lnum_t start_id = tria_in_edges_idx[e_id];
    const cs_lnum_t end_id = tria_in_edges_idx[e_id+1];
    const cs_lnum_t n_tria_in_edge = end_id - start_id;

    if (n_tria_in_edge < 1)
      continue;

    const cs_lnum_t edge_id = selected_edges[e_id];
    const cs_lnum_t v0 = edge_v0[edge_id];
    const cs_lnum_t v1 = v2v_ids[edge_id];
    int n_intx_for_this_edge = 0;

    // edge v0-v1
    cs_real_t e_dir[3] = {vtx_coord[v1][0] - vtx_coord[v0][0],
                          vtx_coord[v1][1] - vtx_coord[v0][1],
                          vtx_coord[v1][2] - vtx_coord[v0][2]};
    cs_real_t e_dir_u[3];
    cs_math_3_normalize(e_dir, e_dir_u);

    // Best candidate of each type
    cs_real_t best_w_up = -1.0, best_w_down = -1.0;
    double best_t_up = 0.0,  best_t_down = 0.0;
    int n_up = 0;
    int n_down = 0;

    /* Loop on triangles intersecting the current edge */
    for (cs_lnum_t tri = start_id; tri < end_id; tri ++) {
      cs_lnum_t tria_id = tria_in_edges_lst[tri];

      const cs_lnum_t ids[3] = {3*tria_id + 0,
                                3*tria_id + 1,
                                3*tria_id + 2};

      const cs_lnum_t vtx_ids_l[3] = {tria_vtx_ids[tria_id][0],
                                      tria_vtx_ids[tria_id][1],
                                      tria_vtx_ids[tria_id][2]};

      cs_real_t coords[3][3];
      for (int dir = 0; dir < 3; dir ++) {
        coords[0][dir] = stl_vtx_coord[vtx_ids_l[0]][dir];
        coords[1][dir] = stl_vtx_coord[vtx_ids_l[1]][dir];
        coords[2][dir] = stl_vtx_coord[vtx_ids_l[2]][dir];
      }

      cs_real_t stl_face_cog_l[3];
      for (cs_lnum_t j = 0; j < 3; j++) {
        // Center of gravity of the STL triangle
        cs_real_t cog_third = (coords[0][j] + coords[1][j] + coords[2][j]);
        stl_face_cog_l[j] = cs_math_1ov3 * cog_third;
      }

      int n_inout[2] = {0, 0};

      /* Compute intersection between the STL plane and the edge v0-v1
         TODO: Not give STL_cog as it is already a triangle
         (not need to cut STL into many triangles as
          it is already a triangle by definition) */
      double t = cs_geom_segment_intersect_face
        (0,
         3,              // n_vtx of the triangle STL
         vtx_ids_l,      // vtx_ids_l, // index of vtx_coord
         stl_vtx_coord,  // stl_vtx_coord, // vtx coordinates of the STL face
         stl_face_cog_l, // center of gravity of the STL face
         vtx_coord[v0],  // begin of the edge v0 < v1
         vtx_coord[v1],  // end of the edge
         n_inout,        // up or down ?
         nullptr);       // Optional : face STL normale

      // Intersection found
      if (t >= 0 && t <= 1) {
        n_intx_for_this_edge += 1;
        intx[edge_id] = t;

        cs_real_t n_unit[3];
        cs_real_t n[3] = {stl_normals[ids[0]],
                          stl_normals[ids[1]],
                          stl_normals[ids[2]]};
        cs_math_3_normalize(n, n_unit);

        /* Reliability criteria "less coplanar" : w = |e.n|.
           If many STL triangles intersects the edge,
           we keep only the intersection with the greatest reliability.
           The intersection can be up or down as explain bellow. */
        cs_real_t w = cs::abs(cs_math_3_dot_product(e_dir_u, n_unit));

        /* Classification from n_inout
           n_inout[0] == 1 && n_inout[1] == 0  => “up”  (+1)
           n_inout[0] == 0 && n_inout[1] == 1  => “down”(-1)
           example: edge cutted : --/-\---/--. In this case, we have respectively
           one up, then down, then up. We retain intersection with the more occurence,
           so two up > one down --> we keep up with the better reliability. */
        int evt = 0;
        if (n_inout[0] == 1 && n_inout[1] == 0) {
          evt = +1;  // up
          n_up += 1; // number of up
        }
        else if (n_inout[0] == 0 && n_inout[1] == 1) {
          evt = -1; // down
          n_down += 1; // number of down
        }
        else {
          assert(0);
        }

        if (evt > 0) {
          if (w > best_w_up) { // less coplanar criteria
            best_w_up = w;
            best_t_up = t;
          }
        }
        else {
          if (w > best_w_down) { // less coplanar criteria
            best_w_down = w;
            best_t_down = t;
          }
        }
      }
    } /* End loop on STL triangle */

    // Go to the next edge if no intersection found
    // or if entering and exiting the edge
    if (   n_intx_for_this_edge == 0
        || (best_w_up < 0.0 && best_w_down < 0.0)
        || n_intx_for_this_edge%2 == 0) {
      intx[edge_id] = -1.0;
      continue;
    }

    // Choose the winner with the best reliability w
    cs_real_t t_best = (n_up > n_down) ? best_t_up : best_t_down;
    cs_real_t w_best = (n_up > n_down) ? best_w_up : best_w_down;

    int s1 = (n_up > n_down) ? VTX_FLUID  : VTX_SOLID;
    int s0 = (s1 == VTX_FLUID) ? VTX_SOLID : VTX_FLUID;

    intx[edge_id] = t_best;
    e_v_idx[edge_id + 1] = 1;

    // Criteria “less coplanar” for each node
    if (w_best > best_angle[v0]) {
      best_angle[v0] = w_best;
      vertex_sign[v0] = s0;
    }
    if (w_best > best_angle[v1]) {
      best_angle[v1] = w_best;
      vertex_sign[v1] = s1;
    }
  }

  CS_FREE(selected_edges);
  CS_FREE(tria_in_edges_idx);
  CS_FREE(tria_in_edges_lst);
  CS_FREE(stl_normals);
  CS_FREE(best_angle);

  /* Add new vertices in the mesh */

  cs_lnum_t n_base_vertices = mesh->n_vertices;
  cs_lnum_t n_add_vtx = 0;

  /* Parallel synchronization */

  cs_gnum_t n_g_edges = 0;
  cs_gnum_t *g_edges_num = nullptr;
  if (cs_glob_n_ranks > 1)
    CS_MALLOC(g_edges_num, n_edges, cs_gnum_t);

  if (cs_glob_n_ranks > 1) {
    n_g_edges = cs::mesh::sync_edges_flag(mesh, v2v,
                                          e_v_idx+1, g_edges_num);
  }
  else
    n_g_edges = n_edges;

  /* Transform counts to index */

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    e_v_idx[i+1] += e_v_idx[i];
    e_v_idx[i] += n_base_vertices;
  }
  e_v_idx[n_edges] += n_base_vertices;

  n_add_vtx = e_v_idx[n_edges] - e_v_idx[0];
  n_base_vertices += n_add_vtx;

  /* Update the vertex sign for new vertices */
  CS_REALLOC(vertex_sign, n_base_vertices, int);
  for (cs_lnum_t i = n_vtx; i < n_vtx + n_add_vtx; i++) {
    vertex_sign[i] = VTX_INTERFACE;
  }

  cs_lnum_t n_vtx_new = mesh->n_vertices + n_add_vtx;
  CS_REALLOC(mesh->vtx_coord, n_vtx_new*3, cs_real_t);
  if (mesh->global_vtx_num != nullptr)
    CS_REALLOC(mesh->global_vtx_num, n_vtx_new, cs_gnum_t);

  _build_edge_vertices(mesh, v2v, intx, n_add_vtx, e_v_idx, g_edges_num);
  cs::mesh::build_add_vertices_gnum(mesh,
                                    n_edges,
                                    n_g_edges,
                                    e_v_idx,
                                    g_edges_num);

  bft_printf(" New vertices of the IBM interface inserted on the main mesh\n");
  bft_printf("   Number of new vertices added: %d\n", n_add_vtx);
  bft_printf("   New total number of vertices: %d\n", n_vtx_new);

  /* Compute n_cut_cell */

  const cs_lnum_t *c2e_idx = c2e->idx;
  const cs_lnum_t *c2e_ids = c2e->ids;

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    const cs_lnum_t s_id = c2e_idx[c_id];
    const cs_lnum_t e_id = c2e_idx[c_id+1];

    int n_intx = 0;

    for (cs_lnum_t eidx = s_id; eidx < e_id; eidx++) {
      cs_lnum_t edge_id = c2e_ids[eidx];
      cs_real_t _intx = intx[edge_id];

      if (_intx > -1) {
        n_intx++;
      }
    }

    if (n_intx >= 3) {
      c_cut[n_c_cut++] = c_id;
    }
  }

  /* Propagates vertex sign only in cutted cells */

  for (cs_lnum_t c_id_loc = 0; c_id_loc < n_c_cut; c_id_loc++) {
    cs_lnum_t c_id = c_cut[c_id_loc];

    bool changed = false;
    bool inconsistency = false;

    const cs_lnum_t s_id = c2e_idx[c_id];
    const cs_lnum_t e_id = c2e_idx[c_id+1];
    const cs_lnum_t n_ec = e_id - s_id;

    cs_lnum_t modified_vtx[24];
    int n_modified_vtx = 0;
    assert(24 >= 2.*n_ec);

    do {
      changed = false;

      for (cs_lnum_t eidx = s_id; eidx < e_id; eidx++) {

        cs_lnum_t edge_id = c2e_ids[eidx];

        if (intx[edge_id] >= 0 && intx[edge_id] <= 1)
          continue;

        cs_lnum_t v0 = edge_v0[edge_id];
        cs_lnum_t v1 = v2v_ids[edge_id];

        int s0 = vertex_sign[v0];
        int s1 = vertex_sign[v1];

        if (s0 != VTX_UNKNOWN && s1 == VTX_UNKNOWN) {
          vertex_sign[v1] = s0;
          changed = true;
          modified_vtx[n_modified_vtx] = v1;
          n_modified_vtx++;
        }
        else if (s1 != VTX_UNKNOWN && s0 == VTX_UNKNOWN) {
          vertex_sign[v0] = s1;
          changed = true;
          modified_vtx[n_modified_vtx] = v0;
          n_modified_vtx++;
        }
        else if (s0 != VTX_UNKNOWN && s1 != VTX_UNKNOWN && s0 != s1) {
#if _DEBUG_
          cs_log_printf
            (CS_LOG_DEFAULT,
             _("Inconsistency in vertex sign: "
               "sign_v0 = %d and sign_v1 = %d should be the same. "
               "Here, a solid region is connected to a fluid region, "
               "whereas an STL surface should delimit the two regions. "
               "This could indicate:\n"
               "- a hole in the STL mesh\n"
               "- the STL mesh is not included in the domain, "
               "so no fluid–solid interface can exist on the boundaries."),
             s0, s1);
#endif
          inconsistency = true;

          /* Disable the cell with inconsistency sign */
          cell_flag[c_id] = 1;

          /* Restore vertex sign of the bad cell */
          for (cs_lnum_t i = 0; i < n_modified_vtx; i++)
            vertex_sign[modified_vtx[i]] = VTX_UNKNOWN;

          break;

        }
      }
    } while (changed && !(inconsistency));
  } /* End loop on cutted cells */

  /* Synchronization for the next loop on internal faces */
  // init i_face_cells in the halo for cs_mesh_init_halo
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (mesh->i_face_cells[f_id][j] >= mesh->n_cells) {
        mesh->i_face_cells[f_id][j] = -1;
        break;
      }
    }
  }

  if (  (mesh->n_domains > 1 || mesh->n_init_perio > 0)
      && mesh->halo == nullptr) {

    cs_halo_type_t halo_type = mesh->halo_type;
    cs_mesh_builder_t *mb = (mesh == cs_glob_mesh) ?
                            cs_glob_mesh_builder :
                            nullptr;
    cs_mesh_init_halo(mesh, mb, halo_type, -1, true);
  }

  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, cell_flag);

  /*
   * Synchronize the signs of vertices shared across MPI interfaces.
   *
   * An inconsistent cell reverts the signs it propagated, while the valid
   * neighboring cell provides the correct face classification. The next
   * operation merges contributions from the owning vertex and its copies,
   * ensuring that all ranks cut the shared face using the signs determined
   * from the valid cell's perspective.
   */

  if (   cs_glob_n_ranks > 1
      && mesh->vtx_interfaces == nullptr) {

    if (mesh->global_vtx_num == nullptr) {
      bft_error(__FILE__, __LINE__, 0,
                "Cannot build vertex interfaces: "
                "mesh->global_vtx_num is null.");
    }

    mesh->vtx_interfaces
      = cs_interface_set_create(mesh->n_vertices,
                                nullptr,
                                mesh->global_vtx_num,
                                mesh->periodicity,
                                0,
                                nullptr,
                                nullptr,
                                nullptr);
  }

  if (mesh->vtx_interfaces != nullptr) {
    cs_datatype_t datatype = cs_datatype_from_type<int>();
    cs_interface_set_inclusive_or(mesh->vtx_interfaces,
                                  mesh->n_vertices,
                                  1,
                                  true,
                                  datatype,
                                  vertex_sign);
  }

  return vertex_sign;
}

/*------------------------------------------------------------------------------
 * A |      | D
 *   |      |
 * B |______| C
 *
 * If a polyline is open, find the start and end vertices (e.g. A and D), and
 * return true.
 * Else, return false.
 *----------------------------------------------------------------------------*/

static bool
_get_open_polyline(cs_lnum_t    fe[],
                   cs_lnum_t    e_stride,
                   uint8_t     *occurs,
                   cs_lnum_t   *start,
                   cs_lnum_t   *end,
                   cs_lnum_2_t  edges[])
{
  for (cs_lnum_t i = 0; i < e_stride; i++) {
    if (fe[i] == -1) continue;
    cs_lnum_t e_id = fe[i];
    occurs[edges[e_id][0]]++;
    occurs[edges[e_id][1]]++;
  }

  *start = -1;
  *end = -1;

  for (cs_lnum_t i = 0; i < e_stride; i++) {
    if (fe[i] == -1) continue;
    cs_lnum_t e_id = fe[i];
    cs_lnum_t p = edges[e_id][0], q = edges[e_id][1];

    if (occurs[p] == 1) {
      if (*start == -1) *start = p; else if (*end == -1) *end = p;
    }
    if (occurs[q] == 1) {
      if (*start == -1) *start = q; else if (*end == -1) *end = q;
    }
  }

  return *start != -1;
}

static inline void
_get_b_face_vertices(cs_mesh_t *mesh,
                     cs_lnum_t  f_id,
                     cs_lnum_t *fv[],
                     cs_lnum_t *nv)
{
  *nv = mesh->b_face_vtx_idx[f_id+1] - mesh->b_face_vtx_idx[f_id];
  *fv = mesh->b_face_vtx_lst + mesh->b_face_vtx_idx[f_id];
}

/*------------------------------------------------------------------------------
 * Remove edge e_id from the face-edge connectivity of face f_id.
 *----------------------------------------------------------------------------*/

static void
_remove_edge_from_face(cs_lnum_t f2e[],
                       cs_lnum_t e_stride,
                       cs_lnum_t f_id,
                       cs_lnum_t e_id,
                       cs_lnum_t e2f[])
{
  cs_lnum_t *fe = f2e + e_stride*f_id;
  cs_lnum_t i;
  for (i = 0; i < e_stride; i++) {
    if (fe[i] == e_id) {
      fe[i] = -1;
      e2f[e_id] = -1;
      break;
    }
  }
  assert(i != e_stride);
}

/*------------------------------------------------------------------------------
 * Add edge e_id to the face-edge connectivity of face f_id.
 *----------------------------------------------------------------------------*/

static void
_add_edge_to_face(cs_lnum_t f2e[],
                  cs_lnum_t e_stride,
                  cs_lnum_t f_id,
                  cs_lnum_t e_id,
                  cs_lnum_t e2f[])
{
  cs_lnum_t *fe = f2e + e_stride*f_id;
  cs_lnum_t i;
  for (i = 0; i < e_stride; i++) {
    if (fe[i] == -1) {
      fe[i] = e_id;
      e2f[e_id] = f_id;
      break;
    }
  }
  assert(i != e_stride);
}

/*------------------------------------------------------------------------------
 * Lexicographical comparison of two vertex coordinates.
 *----------------------------------------------------------------------------*/

static int
_cmp_crd(const cs_real_t x[3], const cs_real_t y[3])
{
  for (int i = 0; i < 3; i++) {
    if (x[i] < y[i]) return -1;
    if (x[i] > y[i]) return 1;
  }
  return 0;
}

/*-----------------------------------------------------------------------------
 * Deduce the face-vertex connectivity from the face-edge connectivity.
 *----------------------------------------------------------------------------*/

static void
_fill_connectivity(cs_lnum_t   b_face_vtx_lst[],
                   cs_lnum_t   e_stride,
                   cs_lnum_t   f_id,
                   cs_lnum_t   fe[],
                   cs_lnum_2_t edges[],
                   cs_lnum_t   b_face_vtx_idx[],
                   cs_lnum_t   ne)
{
  cs_lnum_t *fv = b_face_vtx_lst + e_stride*f_id;
  for (cs_lnum_t j = 0; j < e_stride; j++) fv[j] = -1;

  cs_lnum_t count;
  for (count = 0; count < 2; count++) fv[count] = edges[fe[0]][count];

  while (count != ne) {
    for (cs_lnum_t j = 1; j < e_stride; j++) {
      cs_lnum_t e_id = fe[j];
      if (fe[j] == -1) continue;

      for (int k = 0; k < 2; k++) {
        if (edges[e_id][k] == fv[count-1]) {
          fv[count++] = edges[e_id][(k+1)%2];
          fe[j] = -1;
          break;
        }
      }
      if (fe[j] == -1) break;
    }
  }

  b_face_vtx_idx[f_id+1] = ne;
}

/*------------------------------------------------------------------------------
 * Compute the normal of a face given its vertex connectivity.
 * For now, approximates it with the normal of the triangle formed by the
 * first three vertices.
 * TODO: more robust computation.
 *----------------------------------------------------------------------------*/

static void
_get_normal(const cs_lnum_t  *fv,
            const cs_real_t  *crd,
            cs_real_t         normal[3])
{
  /* Get the normal to the first triangle. */
  cs_lnum_t a = fv[0], b = fv[1], c = fv[2];
  cs_real_t u[3], v[3];
  for (int i = 0; i < 3; i++) {
    u[i] = crd[3*b+i] - crd[3*a+i];
    v[i] = crd[3*c+i] - crd[3*a+i];
  }
  cs_math_3_cross_product(u, v, normal);
}

/*-----------------------------------------------------------------------------
 * After getting the face-vertex connectivity of a new child face,
 * recover its correct orientation using the orientation of its parent face.
 *----------------------------------------------------------------------------*/

static void
_fill_connectivity_and_reorient(cs_mesh_t   *mesh,
                                cs_lnum_t    b_face_vtx_lst[],
                                cs_lnum_t    e_stride,
                                cs_lnum_t    f_id,
                                cs_lnum_t    fe[],
                                cs_lnum_2_t  edges[],
                                cs_lnum_t    b_face_vtx_idx[],
                                cs_lnum_t    ne,
                                cs_lnum_t    o_fid)
{
  _fill_connectivity(b_face_vtx_lst,
                     e_stride,
                     f_id,
                     fe,
                     edges,
                     b_face_vtx_idx,
                     ne);

  /* Get the original normal. */
  cs_real_t o_normal[3] = {0, 0, 0};
  cs_lnum_t *fv = mesh->b_face_vtx_lst + mesh->b_face_vtx_idx[o_fid];
  _get_normal(fv, mesh->vtx_coord, o_normal);

  /* Compute current normal. */
  cs_real_t c_normal[3] = {0.0, 0.0, 0.0};
  fv = b_face_vtx_lst + f_id*e_stride;
  _get_normal(fv, mesh->vtx_coord, c_normal);

  /* Swap the orientation of the current face if necessary. */
  cs_real_t dp = cs_math_3_dot_product(o_normal, c_normal);
  if (dp >= 0.0) return;
  _reverse_array(fv+1, ne-1);
}

/*-----------------------------------------------------------------------------
 * Allocate the buffers necessary for the algorithm with enough memory to be
 * reused across all the cells.
 *----------------------------------------------------------------------------*/

static _cut_data
_allocate_scratch_data(cs_mesh_t *mesh,
                       cs_lnum_t  max_nf,
                       cs_lnum_t  max_nv,
                       cs_lnum_t  n_new_vertices,
                       cs_lnum_t  n_new_cells)
{
  _cut_data cd;
  CS_MALLOC(cd.sd, mesh->n_vertices + n_new_vertices, cs_real_t);
  CS_MALLOC(cd.occurs, mesh->n_vertices + n_new_vertices, uint8_t);
  cd.e_stride = CS_MAX(max_nv+1, max_nf);
  cd.f_size = (2*max_nf+2)*cd.e_stride;
  CS_MALLOC(cd.edges, cd.f_size, cs_lnum_2_t);
  CS_MALLOC(cd.e2f, cd.f_size, cs_lnum_t);
  CS_MALLOC(cd.f2e, cd.f_size, cs_lnum_t);
  CS_MALLOC(cd.xyz, 3*n_new_vertices, cs_real_t);
  CS_MALLOC(cd.indices, n_new_vertices, cs_lnum_t);
  CS_MALLOC(cd.compact, n_new_vertices, cs_lnum_t);
  cd.n_polys = 0;
  CS_MALLOC(cd.polys, 2*n_new_cells, cs_lnum_t);
  return cd;
}

/*-----------------------------------------------------------------------------
 * Cut the cell c_id with the plane (p_normal, p_origin).
 * A cut cell is replaced by the cell lying in the negative half-space
 * spanned by the plane. The cell lying in the positive half-space is appended
 * to the mesh connectivity arrays.
 * Set the new face-vertex connectivity data inside b_face_vtx_idx and
 * b_face_vtx_ids.
 *----------------------------------------------------------------------------*/

static void
_cut_cell(cs_mesh_t            *mesh,
          _cut_data            *cd,
          cs_lnum_t             c_id,
          const cs_real_t       p_normal[],
          const cs_real_t       p_origin[],
          const cs_adjacency_t *c2v,
          const cs_adjacency_t *c2f_b,
          cs_lnum_t            *b_face_vtx_idx,
          cs_lnum_t            *b_face_vtx_lst)
{
  cs_real_t *sd = cd->sd;
  uint8_t *occurs = cd->occurs;
  cs_lnum_t e_stride = cd->e_stride;
  cs_lnum_t f_size = cd->f_size;
  cs_lnum_2_t *edges = cd->edges;
  cs_lnum_t *e2f = cd->e2f;
  cs_lnum_t *f2e = cd->f2e;
  cs_real_t *xyz = cd->xyz;
  cs_lnum_t *indices = cd->indices;
  cs_lnum_t *compact = cd->compact;
  cs_lnum_t *polys = cd->polys;

  /* Reset. */
  for (cs_lnum_t i = 0; i < f_size; i++) {
    edges[i][0] = edges[i][1] = -1;
    e2f[i] = -1;
    f2e[i] = -1;
  }

  /* Compute the signed distances of the cell vertices to the cut plane. */
  int positive = 0, negative = 0;

  for (cs_lnum_t i = c2v->idx[c_id]; i < c2v->idx[c_id+1]; i++) {
    cs_lnum_t p = c2v->ids[i];
    const cs_real_t *crd = mesh->vtx_coord + 3*p;
    sd[p] = 0.0;
    for (int j = 0; j < 3; j++) sd[p] += (crd[j] - p_origin[j]) * p_normal[j];
    if (sd[p] >= _plane_tol) positive++;
    else if (sd[p] <= -_plane_tol) negative++;
    else sd[p] = 0.0;
  }

  /* If the plane does not cut the cell, do nothing. */
  if (positive == 0 || negative == 0) return;

  /* Make the edge-face connectivity. */
  cs_lnum_t n_edges = 0;
  cs_lnum_t nf = c2f_b->idx[c_id+1] - c2f_b->idx[c_id];
  const cs_lnum_t *cf = c2f_b->ids + c2f_b->idx[c_id];

  for (cs_lnum_t i = 0; i < nf; i++) {
    cs_lnum_t f_id = cf[i];
    cs_lnum_t *fv, nv;
    _get_b_face_vertices(mesh, f_id, &fv, &nv);

    cs_lnum_t *fe = f2e + e_stride*i;

    for (cs_lnum_t j = 0; j < nv; j++) {
      cs_lnum_t p = fv[j];
      cs_lnum_t q = fv[(j+1)%nv];
      edges[n_edges][0] = p;
      edges[n_edges][1] = q;
      e2f[n_edges] = i;
      fe[j] = n_edges;
      n_edges++;
    }
  }

  /* Process the edges. */
  cs_lnum_t n_cut_e = 0;

  for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
    cs_lnum_t p = edges[e_id][0];
    cs_lnum_t q = edges[e_id][1];

    cs_real_t dp = sd[p];
    cs_real_t dq = sd[q];

    cs_lnum_t f = e2f[e_id];

    if (dp <= 0 && dq <= 0) {
      /* Edge on the negative side. Do nothing. */
      continue;
    }

    if (dp >= 0 && dq >= 0) {
      /* Edge is on the positive side. */
      _remove_edge_from_face(f2e, e_stride, f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, f + nf, e_id, e2f);
      continue;
    }

    /* Split the edge by the plane. */
    cs_lnum_t a = p, b = q;
    if (a > b) { cs_lnum_t t = a; a = b; b = t; }
    cs_real_t da = sd[a];
    cs_real_t db = sd[b];
    cs_real_t t = da / (da - db);
    cs_real_t *crd = xyz + 3*n_cut_e;
    const cs_real_t *a_crd = mesh->vtx_coord + 3*a;
    const cs_real_t *b_crd = mesh->vtx_coord + 3*b;
    for (int k = 0; k < 3; k++) crd[k] = (1.0-t)*a_crd[k] + t*b_crd[k];

    cs_lnum_t new_e_id = n_edges + n_cut_e;

    if (dp > 0) {
      edges[e_id][0] = q;
      edges[e_id][1] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][0] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][1] = p;
    }
    else {
      edges[e_id][0] = p;
      edges[e_id][1] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][0] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][1] = q;
    }

    _add_edge_to_face(f2e, e_stride, f + nf, new_e_id, e2f);

    n_cut_e++;
  }

  n_edges += n_cut_e;

  /* Get rid of duplicate vertices */
  for (cs_lnum_t i = 0; i < n_cut_e; i++) indices[i] = i;

  cs_lnum_t n_dups = 0;
  for (cs_lnum_t i = 0; i < n_cut_e; i++) {
    cs_real_t *i_crd = xyz + 3*i;
    for (cs_lnum_t j = i+1; j < n_cut_e; j++) {
      cs_real_t *j_crd = xyz + 3*j;
      if (_cmp_crd(i_crd, j_crd) == 0) {
        indices[j] = indices[i];
        n_dups++;
      }
    }
  }

  cs_lnum_t unique_vtx = 0;
  for (int i = 0; i < n_cut_e; i++) compact[i] = -1;

  for (cs_lnum_t i = 0; i < n_cut_e; i++) {
    if (indices[i] == i) {
      compact[i] = unique_vtx++;

      /* Insert the coordinates into the mesh. */
      cs_real_t *crd = mesh->vtx_coord + 3*(mesh->n_vertices+compact[i]);
      const cs_real_t *_crd = xyz + 3*i;
      for (int k = 0; k < 3; k++) crd[k] = _crd[k];
    }
  }

  assert(unique_vtx == n_cut_e-n_dups);

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    for (int j = 0; j < 2; j++) {
      cs_lnum_t p = edges[i][j];
      if (p >= mesh->n_vertices) {
        cs_lnum_t delta = p - mesh->n_vertices;
        cs_lnum_t first = indices[delta];
        cs_lnum_t index = compact[first];
        assert(index != -1);
        edges[i][j] = mesh->n_vertices + index;
      }
    }
  }

  /* Process the faces. */

  cs_lnum_t n_close_e = 0;
  cs_lnum_t seeds[2] = {-1, -1};

  for (cs_lnum_t i = 0; i < nf; i++) {
    cs_lnum_t *fe = f2e + i*e_stride;

    cs_lnum_t old_f = i;
    cs_lnum_t new_f = i + nf;
    cs_lnum_t closing_f = 2*nf;

    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] == -1) continue;
      cs_lnum_t e_id = fe[j];
      occurs[edges[e_id][0]] = 0;
      occurs[edges[e_id][1]] = 0;
    }

    cs_lnum_t start, end;
    if (_get_open_polyline(fe, e_stride, occurs, &start, &end, edges)) {
      /* Polyline is open, close it. */
      cs_lnum_t e_id = n_edges + n_close_e;
      edges[e_id][0] = start;
      edges[e_id][1] = end;

      /* Add the edge to the old and new faces and to the two closing faces. */
      _add_edge_to_face(f2e, e_stride, old_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, new_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, closing_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, closing_f+1, e_id, e2f);

      n_close_e++;

      if (seeds[0] == -1) {
        assert(seeds[1] == -1);
        /* Link the closing polygons with the first faces that share one of its
         * edges. Useful to set their correct orientations later.
         */
        seeds[0] = old_f;
        seeds[1] = new_f;
      }
    }
  }

  assert(seeds[0] != -1);
  assert(seeds[1] != -1);

  /* Fill in the new face connectivities */

  cs_lnum_t f_incr = 0;
  bool seeds_updated = false;

  for (cs_lnum_t i = 0; i < nf; i++) {
    /* Old face. */
    cs_lnum_t *fe = f2e + i*e_stride;

    cs_lnum_t ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    cs_lnum_t f_id = cf[i];

    if (ne == 0) {
      /* Old face belongs to the positive side.
       * Copy its connectivity but do not increment n_b_faces.
       */
      fe = f2e + (i+nf)*e_stride;
      ne = 0;
      for (cs_lnum_t j = 0; j < e_stride; j++) {
        if (fe[j] != -1) ne++;
      }
      assert(ne != 0);

      _fill_connectivity_and_reorient(mesh,
                                      b_face_vtx_lst,
                                      e_stride,
                                      f_id,
                                      fe,
                                      edges,
                                      b_face_vtx_idx,
                                      ne,
                                      f_id);

      mesh->b_face_cells[f_id] = mesh->n_cells;

      continue;
    }

    assert(ne <= e_stride);

    _fill_connectivity_and_reorient(mesh,
                                    b_face_vtx_lst,
                                    e_stride,
                                    f_id,
                                    fe,
                                    edges,
                                    b_face_vtx_idx,
                                    ne,
                                    f_id);

    /* New face. */
    fe = f2e + (i+nf)*e_stride;

    ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    if (ne == 0) continue;

    cs_lnum_t f_id_new = mesh->n_b_faces + f_incr++;

    if (!seeds_updated) {
      if (i == seeds[0]) {
        seeds_updated = true;
        seeds[0] = f_id;
        seeds[1] = f_id_new;
      }
    }

    _fill_connectivity_and_reorient(mesh,
                                    b_face_vtx_lst,
                                    e_stride,
                                    f_id_new,
                                    fe,
                                    edges,
                                    b_face_vtx_idx,
                                    ne,
                                    f_id);

    assert(mesh->b_face_family[f_id_new] == _default_family_id);
    assert(mesh->b_face_cells[f_id_new] == -1);

    mesh->b_face_cells[f_id_new] = mesh->n_cells;
    mesh->b_face_family[f_id_new] = mesh->b_face_family[f_id];
    if (mesh->have_r_gen)
      mesh->b_face_r_c_idx[f_id_new] = mesh->b_face_r_c_idx[f_id];
  }

  assert(seeds_updated);

  for (cs_lnum_t i = 0; i < 2; i++) {
    cs_lnum_t *fe = f2e + (i+2*nf)*e_stride;

    cs_lnum_t ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    assert(ne != 0);
    assert(ne <= e_stride);

    cs_lnum_t f_id = mesh->n_b_faces + f_incr++;

    _fill_connectivity(b_face_vtx_lst,
                       e_stride,
                       f_id,
                       fe,
                       edges,
                       b_face_vtx_idx,
                       ne);

    if (i == 0) {
      mesh->b_face_cells[f_id] = c_id;
    } else {
      mesh->b_face_cells[f_id] = mesh->n_cells;
    }

    mesh->b_face_family[f_id] = _default_family_id;
    if (mesh->have_r_gen) {
      mesh->b_face_r_c_idx[f_id] = 127;  // Mark for later update.
    }

    polys[cd->n_polys++] = f_id;

    /* Set its correct orientation. */
    bool reverse = false;
    cs_lnum_t *fv = b_face_vtx_lst + f_id * e_stride;

    const cs_lnum_t *_fv = b_face_vtx_lst + seeds[i] * e_stride;
    cs_lnum_t _ne = b_face_vtx_idx[seeds[i]+1];

    bool found = false;
    for (cs_lnum_t j = 0; j < _ne; j++) {
      cs_lnum_t _p = _fv[j];
      cs_lnum_t _q = _fv[(j+1)%_ne];
      for (cs_lnum_t k = 0; k < ne && !found; k++) {
        cs_lnum_t p = fv[k];
        cs_lnum_t q = fv[(k+1)%ne];

        if (p == _p && q == _q) {
          found = true;
          reverse = true;
          break;
        }
        else if (p == _q && q == _p) {
          found = true;
          reverse = false;
          break;
        }
      }
    }
    assert(found);

    if (reverse) {
      _reverse_array(fv+1, ne-1);
    }
  }

  /* Update cell family. */
  assert(mesh->cell_family[mesh->n_cells] == _default_family_id);
  mesh->cell_family[mesh->n_cells] = mesh->cell_family[c_id];

  /* Increment mesh counts. */
  mesh->n_cells++;
  mesh->n_b_cells++;
  mesh->n_b_faces += f_incr;
  mesh->n_vertices += unique_vtx;
}

/*-----------------------------------------------------------------------------
 * Adds the cells lying in the positive/negative half-spaces to a group.
 * Useful for further mesh processing.
 *----------------------------------------------------------------------------*/

static void
_mark_cut_cells(cs_mesh_t *mesh, cs_lnum_t n_new_cells, const cs_lnum_t cells[])
{
  cs_mesh_group_cells_set(mesh, "auto:negative_cells", n_new_cells, cells);

  cs_lnum_t *sel_cells = nullptr;
  CS_MALLOC(sel_cells, n_new_cells, cs_lnum_t);
  cs_lnum_t n_c_ini = mesh->n_cells - n_new_cells;
  for (cs_lnum_t i = 0; i < n_new_cells; i++) {
    sel_cells[i] = n_c_ini + i;
  }
  cs_mesh_group_cells_set(mesh, "auto:positive_cells", n_new_cells, sel_cells);
  CS_FREE(sel_cells);
}

/*-----------------------------------------------------------------------------
 *
 * \brief Prepare the cutting of cell faces.
 *
 * - Determine which cells are to be cut,
 * - Transform interior faces to boundary faces for those cells.
 *
 * \param[in, out]   mesh        pointer to mesh structure.
 * \param[in]        p_normals   new number of local vertices
 * \param[in, out]   n_c_cut     new number of local boundary faces
 * \param[out]       c_cut       new number of local cells
 *----------------------------------------------------------------------------*/

static void
_prepare_cut_cell_faces(cs_mesh_t        *mesh,
                        const cs_real_t   p_normals[][3],
                        cs_lnum_t        &n_c_cut,
                        cs_lnum_t         c_cut[])
{
  cs_real_t *p_norms = nullptr;
  CS_MALLOC(p_norms, mesh->n_cells_with_ghosts, cs_real_t);

  n_c_cut = 0;
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    p_norms[c_id] = cs_math_3_square_norm(p_normals[c_id]);
    if (p_norms[c_id] > 0.0)
      c_cut[n_c_cut++] = c_id;
  }

  if (mesh->halo) {
    cs_halo_sync_untyped(mesh->halo,
                         mesh->halo_type,
                         sizeof(cs_real_t),
                         p_norms);
  }

  /* Transform the cut cells internal faces into boundary faces */
  cs_lnum_t *sel_faces = nullptr;
  CS_MALLOC(sel_faces, mesh->n_i_faces, cs_lnum_t);
  cs_lnum_t n_sel = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (p_norms[mesh->i_face_cells[f_id][j]] > 0.0) {
        sel_faces[n_sel++] = f_id;
        break;
      }
    }
  }

  cs_mesh_group_i_faces_set(mesh,
                            "auto:transformed_internal_faces",
                            n_sel,
                            sel_faces);

  cs_mesh_boundary_insert_with_shared_vertices(mesh,
                                               n_sel,
                                               sel_faces);

  CS_FREE(sel_faces);

  cs_mesh_free_rebuildable(mesh, true);

  // TODO: move this inside cs_mesh_free_rebuildable
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (mesh->i_face_cells[f_id][j] >= mesh->n_cells) {
        mesh->i_face_cells[f_id][j] = -1;
        break;
      }
    }
  }

  CS_FREE(p_norms);
}

/*------------------------------------------------------------------------------
 * Update the mesh face connectivity post cut.
 *----------------------------------------------------------------------------*/

static void
_update_face_connectivity(cs_mesh_t  *mesh,
                          cs_lnum_t   b_face_vtx_idx[],
                          cs_lnum_t   b_face_vtx_lst_size,
                          cs_lnum_t  *b_face_vtx_lst[])
{
  CS_FREE(mesh->b_face_vtx_idx);

  /* Counts to indices */
  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    b_face_vtx_idx[i+1] += b_face_vtx_idx[i];
  }
  mesh->b_face_vtx_idx = b_face_vtx_idx;
  mesh->b_face_vtx_connect_size = b_face_vtx_idx[mesh->n_b_faces];

  CS_FREE(mesh->b_face_vtx_lst);

  cs_lnum_t *_b_face_vtx_lst = nullptr;
  CS_MALLOC(_b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_lnum_t);
  cs_lnum_t j = 0;

  cs_lnum_t *lst = *b_face_vtx_lst;

  for (cs_lnum_t i = 0; i < b_face_vtx_lst_size; i++) {
    if (lst[i] != -1)
      _b_face_vtx_lst[j++] = lst[i];
  }

  assert(j == mesh->b_face_vtx_connect_size);

  CS_FREE(lst);

  mesh->b_face_vtx_lst = _b_face_vtx_lst;
}

/*-----------------------------------------------------------------------------
 * Allocate the new face-vertex connectivity arrays.
 *
 * \param[in, out]   mesh            pointer to mesh structure.
 * \param[out]       b_face_vtx_idx  boundary face->vertices index
 * \param[out]       b_face_vtx_lst  boundary face->vertices list
 * \param[in]        n_new_faces  number of new faces
 * \param[in]        e_stride     maximum number of edges per face.
 *
 * \return size of allocated b_face_vtx_lst connectivity array
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_init_b_face_connectivity_arrays(cs_mesh_t  *mesh,
                                 cs_lnum_t  *b_face_vtx_idx[],
                                 cs_lnum_t  *b_face_vtx_lst[],
                                 cs_lnum_t   n_new_faces,
                                 cs_lnum_t   e_stride)
{
  CS_MALLOC(*b_face_vtx_idx, mesh->n_b_faces+n_new_faces+1, cs_lnum_t);
  cs_lnum_t *idx = *b_face_vtx_idx;

  for (cs_lnum_t i = 0; i < mesh->n_b_faces+n_new_faces+1; i++) {
    idx[i] = 0;
  }
  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    idx[i+1] = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
  }

  cs_lnum_t b_face_vtx_lst_size = (mesh->n_b_faces+n_new_faces)*e_stride;
  CS_MALLOC(*b_face_vtx_lst, b_face_vtx_lst_size, cs_lnum_t);
  cs_lnum_t *ids = *b_face_vtx_lst;

  for (cs_lnum_t i = 0; i < b_face_vtx_lst_size; i++) ids[i] = -1;

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    cs_lnum_t *_fv = ids + i*e_stride;

    cs_lnum_t start = mesh->b_face_vtx_idx[i];
    cs_lnum_t end = mesh->b_face_vtx_idx[i+1];
    cs_lnum_t nv = end - start;
    const cs_lnum_t *fv = mesh->b_face_vtx_lst + start;

    for (cs_lnum_t j = 0; j < e_stride && j < nv; j++) {
      _fv[j] = fv[j];
    }
  }

  return b_face_vtx_lst_size;
}

/*------------------------------------------------------------------------------
 * Resize the mesh connectivity arrays before mesh modification.
 *----------------------------------------------------------------------------*/

static void
_resize_mesh(cs_mesh_t  *mesh,
             cs_lnum_t   n_new_vertices,
             cs_lnum_t   n_new_faces,
             cs_lnum_t   n_new_cells)
{
  CS_REALLOC(mesh->vtx_coord, 3*(mesh->n_vertices+n_new_vertices), cs_real_t);

  CS_REALLOC(mesh->b_face_cells, mesh->n_b_faces+n_new_faces, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_new_faces; i++)
    mesh->b_face_cells[mesh->n_b_faces+i] = -1;

  CS_REALLOC(mesh->cell_family, mesh->n_cells+n_new_cells, cs_lnum_t);

  for (cs_lnum_t i = mesh->n_cells; i < mesh->n_cells+n_new_cells; i++)
    mesh->cell_family[i] = _default_family_id;

  CS_REALLOC(mesh->b_face_family, mesh->n_b_faces+n_new_faces, cs_lnum_t);

  for (cs_lnum_t i = mesh->n_b_faces; i < mesh->n_b_faces+n_new_faces; i++)
    mesh->b_face_family[i] = _default_family_id;

  CS_REALLOC(mesh->b_face_r_c_idx, mesh->n_b_faces+n_new_faces, char);
  for (cs_lnum_t i = mesh->n_b_faces; i < mesh->n_b_faces+n_new_faces; i++) {
    mesh->b_face_r_c_idx[i] = 127; // Will be updated
  }
}

/*------------------------------------------------------------------------------
 * Free scratch data used for cell cut algorithm.
 *
 * \param[in, out]  cd  cut helper structure whose members are to be freed.
 *----------------------------------------------------------------------------*/

static void
_free_cut_data(_cut_data  *cd)
{
  CS_FREE(cd->sd);
  CS_FREE(cd->occurs);
  CS_FREE(cd->xyz);
  CS_FREE(cd->indices);
  CS_FREE(cd->compact);
  CS_FREE(cd->polys);
  CS_FREE(cd->f2e);
  CS_FREE(cd->e2f);
  CS_FREE(cd->edges);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut cells with planes.
 *
 * Each cell can be cut by a single plane, defined by its normal and an origin
 * (i.e. any point in the plane). Cells whose assigned normals are null
 * vectors are not cut.
 *
 * The polygons created by the cut are added to a new group,
 * "auto:closing_polygons".
 *
 * This function should be followed by applying a joining on the group,
 * "auto:transformed_internal_faces".
 *
 * \param[in, out]  mesh      mesh to cut
 * \param[in]       p_normals array of plane_normals of size mesh->n_cells
 * \param[in]       p_origins array of plane origins of size mesh->n_cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cut(cs_mesh_t       *mesh,
            const cs_real_t  p_normals[][3],
            const cs_real_t  p_origins[][3])
{
  std::chrono::high_resolution_clock::time_point t0
    = std::chrono::high_resolution_clock::now();

  bft_printf("\nStart of cell-plane cut\n\n");

  cs_lnum_t *sel_cells = nullptr;
  CS_MALLOC(sel_cells, mesh->n_cells, cs_lnum_t);
  cs_lnum_t n_new_cells;

  _prepare_cut_cell_faces(mesh, p_normals, n_new_cells, sel_cells);

  /* Create useful cell connectivity arrays. */
  const cs_lnum_t n_c_ini = mesh->n_cells;
  const cs_lnum_t n_v_ini = mesh->n_vertices;
  const cs_lnum_t n_b_ini = mesh->n_b_faces;
  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();
  cs_adjacency_t *c2f_b = cs_mesh_adjacency_c2f_boundary(mesh);

  /* Size-up the problem. */
  cs_lnum_t max_nv = 0, max_nf = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t c_id = mesh->b_face_cells[f_id];
    cs_lnum_t nf = c2f_b->idx[c_id+1] - c2f_b->idx[c_id];
    if (nf > max_nf) max_nf = nf;

    cs_lnum_t nv = mesh->b_face_vtx_idx[f_id+1] - mesh->b_face_vtx_idx[f_id];
    if (nv > max_nv) max_nv = nv;
  }

  cs_lnum_t n_new_faces = n_new_cells*(max_nf + 2);
  cs_lnum_t n_new_vertices = (n_new_faces - 2*n_new_cells)*2;

  _cut_data cd = _allocate_scratch_data(mesh,
                                        max_nf,
                                        max_nv,
                                        n_new_vertices,
                                        n_new_cells);

  /* Init b_face connectivity arrays. */
  cs_lnum_t *b_face_vtx_idx = nullptr;
  cs_lnum_t *b_face_vtx_lst = nullptr;
  cs_lnum_t b_face_vtx_lst_size;
  b_face_vtx_lst_size = _init_b_face_connectivity_arrays(mesh,
                                                         &b_face_vtx_idx,
                                                         &b_face_vtx_lst,
                                                         n_new_faces,
                                                         cd.e_stride);

  /* Resize mesh. */
  _resize_mesh(mesh, n_new_vertices, n_new_faces, n_new_cells);

  std::chrono::high_resolution_clock::time_point t1
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_prepare
    = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

  /* Cut. */
  for (cs_lnum_t i = 0; i < n_new_cells; i++) {
    cs_lnum_t c_id = sel_cells[i];
    const cs_real_t *p_normal = p_normals[c_id];
    const cs_real_t *p_origin = p_origins[c_id];
    _cut_cell(mesh,
              &cd,
              c_id,
              p_normal,
              p_origin,
              c2v,
              c2f_b,
              b_face_vtx_idx,
              b_face_vtx_lst);
  }

  std::chrono::high_resolution_clock::time_point t2
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_cut
    = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

  n_new_cells = mesh->n_cells - n_c_ini;

  /* Add the closing polygons to a group. */
  cs_mesh_group_b_faces_set(mesh,
                            "auto:closing_polygons",
                            cd.n_polys,
                            cd.polys);

  /* Free cut data. */
  _free_cut_data(&cd);
  cs_adjacency_destroy(&c2f_b);
  cs_mesh_adjacencies_finalize();

  /* Set up the end b_face-vertex connectivity arrays. */
  _update_face_connectivity(mesh,
                            b_face_vtx_idx,
                            b_face_vtx_lst_size,
                            &b_face_vtx_lst);

  /* Optional: mark the cut cells for further post-processing. */
  _mark_cut_cells(mesh, n_new_cells, sel_cells);
  CS_FREE(sel_cells);

  /* Update parallel data. */
  n_new_vertices = mesh->n_vertices - n_v_ini;
  n_new_faces = mesh->n_b_faces - n_b_ini;
  _update_parallelism(mesh,
                      true, n_new_vertices,
                      true, n_new_faces,
                      false, 0,
                      n_new_cells);
  cs_mesh_update_auxiliary(mesh);

  /* Tag mesh for repartitionning. */
  mesh->modified |= CS_MESH_MODIFIED;
  mesh->modified |= CS_MESH_MODIFIED_BALANCE;

  bft_printf("\nEnd of cell-plane cut.\n");

  t1 = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_update
    = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t2);

  if (mesh->verbosity > 0) {

    cs_mesh_print_element_counts(mesh,
                                 _("Mesh after cells cut"));

    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\nMesh cells cut:\n\n"
         "  Preparation:                                  %.3g\n"
         "  Cells cut:                                    %.3g\n"
         "  Mesh update:                                  %.3g\n"),
       (double)(te_prepare.count()*1.e-6),
       (double)(te_cut.count()*1.e-6),
       (double)(te_update.count()*1.e-6));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut cells edges with STL triangles.
 *
 * Cells mesh are cutted by triangles according to the STL format. The algorithm
 * computes the intersection between edges mesh and STL triangles to conserves
 * conformity but the results can be warped immersed faces.
 *
 * The polygons created by the cut are added to a new group,
 * "auto:closing_polygons".
 *
 * \param[in]       stl_mesh    STL mesh
 * \param[in, out]  mesh        mesh to cut
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cut_edges_by_stl(const char  *stl_file_name,
                         cs_mesh_t   *mesh)
{
  std::chrono::high_resolution_clock::time_point t0
    = std::chrono::high_resolution_clock::now();

  bft_printf("\n Start intersection between the main mesh"
             " and STL triangles (edge version)\n"
             " ================================================="
             "========================\n\n");

  bool remove_cells = true;
  cs_mesh_quantities_t *mq_old = nullptr;

  if (remove_cells) {
    mq_old = cs_mesh_quantities_create();
    cs_mesh_quantities_compute(mesh, mq_old);
  }

  /* Free mesh data that will be rebuilt, as it would become
     inconsistent once the mesh is modified. */
  cs_mesh_free_rebuildable(mesh, true);

  cs_stl_mesh_t *stl_mesh = cs_stl_mesh_add("my_exemple");
  cs_stl_file_read(stl_mesh, stl_file_name);

  /* Get the time control from the default writer (writer_id = -1) */
  //cs_time_control_t *tc = cs_post_get_time_control(-1);

  bool output_at_start = true;
  bool output_at_end = false;
  int interval_n = -1;
  cs_real_t interval_t = -1.0; //tc->interval_t;

  cs_stl_post_init_writer("STL_OBJECTS",
                          "postprocessing",
                          "Ensight Gold",
                          "",
                          FVM_WRITER_FIXED_MESH, //FVM_WRITER_TRANSIENT_COORDS,
                          output_at_start,
                          output_at_end,
                          interval_n,
                          interval_t);

  cs_stl_post_add_mesh(stl_mesh);

  const cs_lnum_t n_vtx = mesh->n_vertices;

  cs_adjacency_t *v2v = cs_mesh_adjacency_v2v(mesh);
  cs_lnum_t *v2v_idx = v2v->idx;
  cs_lnum_t n_edges = v2v_idx[n_vtx];

  cs_adjacency_t *f2e = cs::mesh::build_f2e_connect(mesh, v2v);
  cs_adjacency_t *c2f = cs_mesh_adjacency_c2f(mesh, 1);
  cs_adjacency_t *c2e = cs_adjacency_compose(n_edges, c2f, f2e);

  cs_lnum_t *edge_v0;
  CS_MALLOC(edge_v0, n_edges, cs_lnum_t);

  double *intx;
  CS_MALLOC(intx, n_edges, double);

  int *vertex_sign;
  CS_MALLOC(vertex_sign, n_vtx, int);

  /* For each edge, e_v_idx will contain the starting index of
     vertices inserted on edges requiring subdivision */

  cs_lnum_t *e_v_idx;
  CS_MALLOC(e_v_idx, n_edges+1, cs_lnum_t);
  e_v_idx[0] = 0;

  cs_lnum_t n_cut_cells = 0;
  cs_lnum_t *sel_cells;
  CS_MALLOC(sel_cells, mesh->n_cells, cs_lnum_t);

  int *cell_flag_inconsistency;
  CS_MALLOC(cell_flag_inconsistency, mesh->n_cells_with_ghosts, int);
  cs_lnum_t *n2o_cells;
  CS_MALLOC(n2o_cells, 2*mesh->n_cells, cs_lnum_t);

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    n2o_cells[c_id] = c_id;
    cell_flag_inconsistency[c_id] = 0;
  }

  vertex_sign = _cut_edges(stl_mesh,
                           mesh,
                           v2v,
                           c2e,
                           edge_v0,
                           intx,
                           vertex_sign,
                           e_v_idx,
                           n_cut_cells,
                           sel_cells,
                           cell_flag_inconsistency);

  bft_printf("\n Intersection between edges and STL triangles"
             " successfully computed\n");
  bft_printf(" --> Starting mesh face subdivision\n\n");

  std::chrono::high_resolution_clock::time_point t1
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_cut_edges
    = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

  cs_lnum_t n_i_faces_old = mesh->n_i_faces;

  cs_lnum_t *b_sub_face_vtx_idx, *b_sub_face_vtx_lst;
  cs_lnum_t n_sub_face_max = 10*mesh->n_b_faces;
  cs_lnum_t n_vtx_max = n_sub_face_max*10;
  CS_MALLOC(b_sub_face_vtx_idx, n_sub_face_max + 1, cs_lnum_t);
  CS_MALLOC(b_sub_face_vtx_lst, n_vtx_max, cs_lnum_t);
  b_sub_face_vtx_idx[0] = 0;

  cs_lnum_t *b_face_o2n_idx, *b_face_o2n_connect_idx;
  CS_MALLOC(b_face_o2n_idx, mesh->n_b_faces + 1, cs_lnum_t);
  CS_MALLOC(b_face_o2n_connect_idx, mesh->n_b_faces + 1, cs_lnum_t);
  b_face_o2n_idx[0] = 0; // sub_face index
  b_face_o2n_connect_idx[0] = 0; // vertex index

  cs_lnum_t *i_face_o2n_idx, *i_face_o2n_connect_idx;
  CS_MALLOC(i_face_o2n_idx, mesh->n_i_faces + 1, cs_lnum_t);
  CS_MALLOC(i_face_o2n_connect_idx, mesh->n_i_faces + 1, cs_lnum_t);
  i_face_o2n_idx[0] = 0; // sub_face index
  i_face_o2n_connect_idx[0] = 0; // vertex index

  cs_lnum_t *i_sub_face_vtx_idx, *i_sub_face_vtx_lst;
  cs_lnum_t n_sub_i_face_max = 10*mesh->n_i_faces;
  cs_lnum_t n_i_face_vtx_max = n_sub_i_face_max*10;
  CS_MALLOC(i_sub_face_vtx_idx, n_sub_i_face_max + 1, cs_lnum_t);
  CS_MALLOC(i_sub_face_vtx_lst, n_i_face_vtx_max, cs_lnum_t);
  i_sub_face_vtx_idx[0] = 0;

  cs_lnum_t *dark_edge_in_b_face_idx, *light_edge_in_b_face_idx;
  CS_MALLOC(dark_edge_in_b_face_idx, mesh->n_b_faces+1, cs_lnum_t);
  CS_MALLOC(light_edge_in_b_face_idx, mesh->n_b_faces+1, cs_lnum_t);

  cs_lnum_2_t *dark_edge_in_b_face, *light_edge_in_b_face;
  CS_MALLOC(dark_edge_in_b_face, 10*mesh->n_b_faces, cs_lnum_2_t);
  CS_MALLOC(light_edge_in_b_face, 10*mesh->n_b_faces, cs_lnum_2_t);

  dark_edge_in_b_face_idx[0] = 0;
  light_edge_in_b_face_idx[0] = 0;

  cs_lnum_t *dark_edge_in_i_face_idx, *light_edge_in_i_face_idx;
  CS_MALLOC(dark_edge_in_i_face_idx, mesh->n_i_faces+1, cs_lnum_t);
  CS_MALLOC(light_edge_in_i_face_idx, mesh->n_i_faces+1, cs_lnum_t);

  cs_lnum_2_t *dark_edge_in_i_face, *light_edge_in_i_face;
  CS_MALLOC(dark_edge_in_i_face, 10*mesh->n_i_faces, cs_lnum_2_t);
  CS_MALLOC(light_edge_in_i_face, 10*mesh->n_i_faces, cs_lnum_2_t);

  dark_edge_in_i_face_idx[0] = 0;
  light_edge_in_i_face_idx[0] = 0;

  int n_sub_i_face_tot = 0, n_i_face_vtx_tot = 0;
  int n_sub_b_face_tot = 0, n_b_face_vtx_tot = 0;

  _cut_faces(mesh,
             intx,
             v2v,
             f2e,
             edge_v0,
             vertex_sign,
             cell_flag_inconsistency,
             e_v_idx,
             n_sub_b_face_tot,
             n_b_face_vtx_tot,
             n_sub_i_face_tot,
             n_i_face_vtx_tot,
             b_face_o2n_idx,
             b_face_o2n_connect_idx,
             i_face_o2n_idx,
             i_face_o2n_connect_idx,
             b_sub_face_vtx_idx,
             b_sub_face_vtx_lst,
             i_sub_face_vtx_idx,
             i_sub_face_vtx_lst,
             light_edge_in_b_face_idx,
             light_edge_in_i_face_idx,
             light_edge_in_b_face,
             light_edge_in_i_face,
             dark_edge_in_b_face_idx,
             dark_edge_in_i_face_idx,
             dark_edge_in_b_face,
             dark_edge_in_i_face);

  cs_lnum_t n_b_faces_tot = b_face_o2n_idx[mesh->n_b_faces];
  cs_lnum_t *_b_face_cells;
  CS_MALLOC(_b_face_cells, n_b_faces_tot, cs_lnum_t);

  cs_lnum_t n_i_faces_tot = i_face_o2n_idx[mesh->n_i_faces];
  cs_lnum_2_t *_i_face_cells;
  CS_MALLOC(_i_face_cells, n_i_faces_tot, cs_lnum_2_t);

  _update_b_face_connectivity(mesh,
                              n_b_face_vtx_tot,
                              b_face_o2n_idx,
                              b_sub_face_vtx_idx,
                              b_sub_face_vtx_lst,
                              _b_face_cells);

  _update_i_face_connectivity(mesh,
                              n_i_face_vtx_tot,
                              i_face_o2n_idx,
                              i_sub_face_vtx_idx,
                              i_sub_face_vtx_lst,
                              _i_face_cells);

  bft_printf(" Face connectivity successfully rebuilt\n");
  bft_printf(" --> Starting cell mesh subdivision\n\n");

  std::chrono::high_resolution_clock::time_point t2
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_cut_faces
    = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

  /* Realloc to add the new immersed plane computed in _cut_cells
     Generally, one plane cuts the cells into two new cells
     Multiply by two for security (for multiple cycle in one cell) */

  const cs_real_t n_ib_plane_max = 2. * n_cut_cells;
  CS_REALLOC(mesh->i_face_vtx_idx, n_i_faces_tot + n_ib_plane_max + 1, cs_lnum_t);
  CS_REALLOC(_i_face_cells, n_i_faces_tot + n_ib_plane_max, cs_lnum_2_t);

  /* 10 points maximum for each immersed plane */

  const cs_real_t n_ib_vtx_max = 10. * n_ib_plane_max;
  CS_REALLOC(mesh->i_face_vtx_lst, n_i_face_vtx_tot + n_ib_vtx_max, cs_lnum_t);
  CS_REALLOC(mesh->cell_family, mesh->n_cells + n_ib_plane_max, cs_lnum_t);

  const cs_lnum_t n_i_faces_tot_without_ibm_plane = n_i_faces_tot;
  cs_lnum_t n_new_cells = 0;

  _cut_cells(mesh,
             n_i_faces_old,
             n_new_cells,
             n_cut_cells,
             sel_cells,
             c2f,
             vertex_sign,
             e_v_idx,
             cell_flag_inconsistency,
             b_face_o2n_idx,
             i_face_o2n_idx,
             b_sub_face_vtx_idx,
             b_sub_face_vtx_lst,
             i_sub_face_vtx_idx,
             i_sub_face_vtx_lst,
             light_edge_in_b_face_idx,
             light_edge_in_i_face_idx,
             light_edge_in_b_face,
             light_edge_in_i_face,
             n2o_cells,
             n_i_faces_tot,
             n_i_face_vtx_tot,
             _b_face_cells,
             _i_face_cells);

  const cs_lnum_t n_i_faces_tot_with_ibm_plane = n_i_faces_tot;
  const cs_lnum_t n_ib_plane
    =  n_i_faces_tot_with_ibm_plane - n_i_faces_tot_without_ibm_plane;

  /* Ajust reallocation */
  CS_REALLOC(mesh->i_face_vtx_idx, n_i_faces_tot + 1, cs_lnum_t);
  CS_REALLOC(_i_face_cells, n_i_faces_tot, cs_lnum_2_t);
  CS_REALLOC(mesh->i_face_vtx_lst, n_i_face_vtx_tot, cs_lnum_t);
  CS_REALLOC(mesh->i_face_r_gen, n_i_faces_tot, char);
  CS_REALLOC(mesh->i_face_family, n_i_faces_tot, int);
  CS_REALLOC(mesh->cell_family, mesh->n_cells + n_new_cells, cs_lnum_t);

  for (cs_lnum_t i = mesh->n_i_faces; i < n_i_faces_tot; i++) {
    mesh->i_face_r_gen[i] = 0;
    mesh->i_face_family[i] = 1;
  }
  mesh->i_face_vtx_connect_size = mesh->i_face_vtx_idx[n_i_faces_tot];

  cs_lnum_t n_cells_old = mesh->n_cells;
  mesh->n_cells += n_new_cells;
  mesh->n_i_faces = n_i_faces_tot;
  mesh->n_b_faces = n_b_faces_tot;

  CS_FREE(mesh->b_face_cells);
  mesh->b_face_cells = _b_face_cells;

  CS_FREE(mesh->i_face_cells);
  mesh->i_face_cells = _i_face_cells;

  _update_parallelism(mesh,
                      false, 0,
                      false, 0,
                      true, n_ib_plane,
                      n_new_cells);

  bft_printf(" Cells successfully subdivided\n");
  bft_printf(" --> Starting transformation of each internal IBM interface"
             "\n     into two boundary faces\n\n");

  std::chrono::high_resolution_clock::time_point t3
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_cut_cells
    = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);

  /* Convert each interior IBM faces into two boundary faces */

  const cs_lnum_t n_i_faces_without_ibm_new = mesh->n_i_faces - n_ib_plane;
  const cs_lnum_t n_b_faces_old = mesh->n_b_faces;

  cs_lnum_t *ib_face_id;
  CS_MALLOC(ib_face_id, n_ib_plane, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_ib_plane; i++) {
    ib_face_id[i] = n_i_faces_without_ibm_new + i;
  }

  cs_mesh_boundary_insert_with_shared_vertices(mesh,
                                               n_ib_plane,
                                               ib_face_id);

  CS_FREE(ib_face_id);

  const cs_lnum_t n_b_faces_new = mesh->n_b_faces - n_b_faces_old;
  //const cs_lnum_t n_vertices_new = mesh->n_vertices - n_vertices_old;

  cs_lnum_t *ib_face_id_new;
  CS_MALLOC(ib_face_id_new, n_b_faces_new, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_b_faces_new; i++) {
    ib_face_id_new[i] = n_b_faces_old + i;
  }

  // Add the closing polygons to a group.
  cs_mesh_group_b_faces_add(mesh,
                            "auto:closing_polygons",
                            n_b_faces_new,
                            ib_face_id_new);

  cs_boundary_zone_define("ibm_faces", "auto:closing_polygons", 0);

  mesh->modified |= (CS_MESH_MODIFIED | CS_MESH_MODIFIED_BALANCE);

  CS_FREE(ib_face_id_new);

  bft_printf(" Internal to boundary transformation step completed\n");

  if (remove_cells) {
    bft_printf(" Start to remove the solid part of the mesh\n");

    _remove_invalid_cells(mesh,
                          mq_old,
                          n_cells_old,
                          v2v,
                          c2e,
                          edge_v0,
                          cell_flag_inconsistency,
                          n2o_cells,
                          vertex_sign);

    cs_mesh_quantities_destroy(mq_old);

    bft_printf(" Solid part of the mesh successfuly removed on the mesh\n");

    std::chrono::high_resolution_clock::time_point t4
      = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds te_remove_cells
      = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t0);
  }

  cs_renumber_b_faces(cs_glob_mesh);
  cs_renumber_i_faces(cs_glob_mesh);

  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);

  cs_mesh_update_auxiliary(cs_glob_mesh);
  cs_gradient_free_quantities();
  cs_matrix_update_mesh();
  cs_mesh_update_selectors(mesh);
  cs_mesh_location_build(mesh, -1);

  /* Destroy */
  cs_adjacency_destroy(&c2e);
  cs_adjacency_destroy(&c2f);
  cs_adjacency_destroy(&f2e);
  cs_adjacency_destroy(&v2v);

  CS_FREE(e_v_idx);
  CS_FREE(vertex_sign);
  CS_FREE(edge_v0);
  CS_FREE(intx);

  CS_FREE(b_face_o2n_idx);
  CS_FREE(i_face_o2n_idx);
  CS_FREE(b_face_o2n_connect_idx);
  CS_FREE(i_face_o2n_connect_idx);
  CS_FREE(b_sub_face_vtx_idx);
  CS_FREE(b_sub_face_vtx_lst);
  CS_FREE(i_sub_face_vtx_idx);
  CS_FREE(i_sub_face_vtx_lst);

  CS_FREE(light_edge_in_b_face_idx);
  CS_FREE(light_edge_in_i_face_idx);
  CS_FREE(light_edge_in_b_face);
  CS_FREE(light_edge_in_i_face);
  CS_FREE(dark_edge_in_b_face_idx);
  CS_FREE(dark_edge_in_i_face_idx);
  CS_FREE(dark_edge_in_b_face);
  CS_FREE(dark_edge_in_i_face);

  CS_FREE(n2o_cells);
  CS_FREE(sel_cells);
  CS_FREE(cell_flag_inconsistency);

  { /* Check the single_faces_to_cells */

    cs_mesh_adjacencies_initialize();
    cs_mesh_adjacencies_update_mesh();
    //cs_mesh_adjacencies_update_cell_i_faces();

    const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
    bft_printf(" single_faces_to_cells = %d\n",
               ma->single_faces_to_cells);
    if (ma->single_faces_to_cells == false)
     cs_log_warning(_("%s: the mesh have not single faces to cells."
                      " The code is not designed to...\n"),
                    __func__);
    cs_mesh_adjacencies_finalize();
  }

  if (mesh->verbosity > 0) {

    cs_mesh_print_element_counts(mesh,
                                 _(" Mesh after cells cut"));

    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);

    /*cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\nMesh cells cut:\n\n"
         "  Preparation:                                  %.3g\n"
         "  Cells cut:                                    %.3g\n"
         "  Mesh update:                                  %.3g\n"),
       (double)(te_prepare.count()*1.e-6),
       (double)(te_cut.count()*1.e-6),
       (double)(te_update.count()*1.e-6));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);*/

  }
}

/*----------------------------------------------------------------------------*/
