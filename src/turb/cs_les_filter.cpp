/*============================================================================
 * Filters for dynamic models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "alge/cs_cell_to_vertex.h"
#include "base/cs_dispatch.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_les_filter.h"

/*----------------------------------------------------------------------------*/

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

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute filters for dynamic models.
 *
 *
 * This function deals with the standard or extended neighborhood.
 *
 * \param[in]   ctx      reference to dispatch context
 * \param[in]   val      array of values to filter
 * \param[out]  f_val    array of filtered values
 */
/*----------------------------------------------------------------------------*/

static void
_les_filter_ext_neighborhood_scalar(cs_dispatch_context  &ctx,
                                    const cs_real_t       val[],
                                    cs_real_t             f_val[])

{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t  *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t  *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2c = ma->cell_cells;
  if (c2c == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
  }

  assert(cell_cells_idx != nullptr);

  /* Synchronize variable */

  cs_halo_sync_var(mesh->halo, CS_HALO_EXTENDED, const_cast<cs_real_t *>(val));

  /* Define filtered variable array */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Initialization */
    cs_real_t w1 = val[c_id] * cell_vol[c_id];
    cs_real_t w2 = cell_vol[c_id];

    /* Contribution from extended neighborhood */

    const cs_lnum_t s_id_ex = cell_cells_idx[c_id];
    const cs_lnum_t e_id_ex = cell_cells_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_ex; cidx < e_id_ex; cidx++) {
      cs_lnum_t c_id_adj = cell_cells_lst[cidx];

      w1 += val[c_id_adj] * cell_vol[c_id_adj];
      w2 += cell_vol[c_id_adj];
    }

    /* Contribution from adjacent faces */

    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      cs_lnum_t c_id_adj = c2c[cidx];

      w1 += val[c_id_adj] * cell_vol[c_id_adj];
      w2 += cell_vol[c_id_adj];
    }

    f_val[c_id] = w1 / w2;

  }); /* End of loop on cells */

  /* Synchronize variable */

  ctx.wait();
  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_val);
}

END_C_DECLS

/*============================================================================
 * Private C++ function definitions
 *============================================================================*/

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute filters for dynamic models.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * This function deals with the standard or extended neighborhood.
 *
 * \param[in]   ctx      reference to dispatch context
 * \param[in]   val      array of values to filter
 * \param[out]  f_val    array of filtered values
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_les_filter_ext_neighborhood_strided(cs_dispatch_context  &ctx,
                                     const cs_real_t       val[][stride],
                                     cs_real_t             f_val[][stride])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t *cell_cells_idx = mesh->cell_cells_idx;
  const cs_lnum_t *cell_cells_lst = mesh->cell_cells_lst;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2c = ma->cell_cells;
  if (c2c == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
  }

  using var_t = cs_real_t[stride];

  assert(cell_cells_idx != nullptr);

  /* Allocate and initialize working buffers */

  /* Synchronize variable */

  cs_halo_sync_r(mesh->halo, CS_HALO_EXTENDED, ctx.use_gpu(),
                 const_cast<var_t *>(val));

  /* Define filtered variable array */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    var_t w1, w2;

    for (cs_lnum_t j = 0; j < stride; j++) {
      w1[j] = val[c_id][j] * cell_vol[c_id];
      w2[j] = cell_vol[c_id];
    }

    /* Contribution from extended neighborhood */

    const cs_lnum_t s_id_ex = cell_cells_idx[c_id];
    const cs_lnum_t e_id_ex = cell_cells_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_ex; cidx < e_id_ex; cidx++) {
      cs_lnum_t c_id_adj = cell_cells_lst[cidx];

      for (cs_lnum_t j = 0; j < stride; j++) {
        w1[j] += val[c_id_adj][j] * cell_vol[c_id_adj];
        w2[j] += cell_vol[c_id_adj];
      }
    }

    /* Contribution from adjacent faces */

    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      cs_lnum_t c_id_adj = c2c[cidx];

      for (cs_lnum_t j = 0; j < stride; j++) {
        w1[j] += val[c_id_adj][j] * cell_vol[c_id_adj];
        w2[j] += cell_vol[c_id_adj];
      }
    }

    for (cs_lnum_t j = 0; j < stride; j++)
      f_val[c_id][j] = w1[j]/w2[j];

  }); /* End of loop on cells */

  /* Synchronize variable */

  ctx.wait();
  cs_halo_sync_r(mesh->halo, CS_HALO_EXTENDED, ctx.use_gpu(), f_val);
}

#endif /* cplusplus */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute filters for dynamic models.
 *
 *
 *
 * This function deals with the standard or extended neighborhood.
 *
 * \param[in]   val      array of values to filter
 * \param[out]  f_val    array of filtered values
 */
/*----------------------------------------------------------------------------*/

void
cs_les_filter_scalar(const cs_real_t  val[],
                     cs_real_t        f_val[])
{
  cs_dispatch_context ctx;

  if (cs_ext_neighborhood_get_type() == CS_EXT_NEIGHBORHOOD_COMPLETE) {
    _les_filter_ext_neighborhood_scalar(ctx, val, f_val);
    return;
  }

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t  n_vertices = mesh->n_vertices;

  /* Allocate and initialize working buffer */

  cs_real_t *v_val, *v_weight;
  CS_MALLOC_HD(v_val, n_vertices, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(v_weight, n_vertices, cs_real_t, cs_alloc_mode);

  /* Define filtered variable array */

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0,
                    1,    /* v_val dimension */
                    true, /* ignore periodicity of rotation */
                    cell_vol,
                    val,
                    nullptr,
                    v_val);

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0,
                    1,
                    true, /* ignore periodicity of rotation */
                    nullptr,
                    cell_vol,
                    nullptr,
                    v_weight);

  /* Build cell average */

  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();
  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    cs_real_t _f_val = 0, _f_weight = 0;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      _f_val += v_val[v_id] * v_weight[v_id];
      _f_weight += v_weight[v_id];
    }

    f_val[c_id] = _f_val / _f_weight;
  });

  ctx.wait();

  CS_FREE_HD(v_weight);
  CS_FREE_HD(v_val);

  /* Synchronize variable */
  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, ctx.use_gpu(), f_val);
}

END_C_DECLS

/*=============================================================================
 * Public C++ function definitions
 *============================================================================*/

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute filters for dynamic models.
 *
 * template parameters:
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * This function deals with the standard or extended neighborhood.
 *
 * \param[in]   val      array of values to filter
 * \param[out]  f_val    array of filtered values
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_les_filter_strided(const cs_real_t  val[][stride],
                      cs_real_t        f_val[][stride])
{
  assert(stride == 3 or stride == 6);

  cs_dispatch_context ctx;

  if (cs_ext_neighborhood_get_type() == CS_EXT_NEIGHBORHOOD_COMPLETE) {
    _les_filter_ext_neighborhood_strided<stride>(ctx, val, f_val);
    return;
  }

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_lnum_t  n_vertices = mesh->n_vertices;

  using var_t = cs_real_t[stride];

  /* Allocate and initialize working buffer */
  var_t *v_val;
  cs_real_t *v_weight;
  CS_MALLOC_HD(v_val, n_vertices, var_t, cs_alloc_mode);
  CS_MALLOC_HD(v_weight, n_vertices, cs_real_t, cs_alloc_mode);

  /* Define filtered variable array */

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0,
                    stride,
                    true, /* ignore periodicity of rotation */
                    cell_vol,
                    reinterpret_cast<const cs_real_t *>(val),
                    nullptr,
                    reinterpret_cast<cs_real_t *>(v_val));

  cs_cell_to_vertex(CS_CELL_TO_VERTEX_LR,
                    0,
                    1,
                    true, /* ignore periodicity of rotation */
                    nullptr,
                    cell_vol,
                    nullptr,
                    v_weight);

  /* Build cell average */

  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();
  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    cs_real_t _f_val[9], _f_weight = 0.;

    for (cs_lnum_t k = 0; k < stride; k++)
      _f_val[k] = 0;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];

      for (cs_lnum_t k = 0; k < stride; k++) {
        _f_val[k] += v_val[v_id][k] * v_weight[v_id];
      }
      _f_weight += v_weight[v_id];
    }

    for (cs_lnum_t k = 0; k < stride; k++) {
      f_val[c_id][k] = _f_val[k] / _f_weight;
    }

  });

  ctx.wait();

  CS_FREE_HD(v_weight);
  CS_FREE_HD(v_val);

  /* Synchronize variable */
  cs_halo_sync_r(mesh->halo, CS_HALO_STANDARD, ctx.use_gpu(), f_val);
}

// Force instanciation

template void
cs_les_filter_strided(const cs_real_t  val[][3],
                      cs_real_t        f_val[][3]);

template void
cs_les_filter_strided(const cs_real_t  val[][6],
                      cs_real_t        f_val[][6]);

/*----------------------------------------------------------------------------*/

#endif /* cplusplus */
