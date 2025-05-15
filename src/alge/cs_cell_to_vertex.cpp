/*============================================================================
 * Cell to vertex interpolation.
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

#include "alge/cs_blas.h"
#include "base/cs_array.h"
#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_cell_to_vertex.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cell_to_vertex.cpp
 * \brief Cell to vertex interpolation..
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.E+12
#endif

/*=============================================================================
 * Local type definitions
 *============================================================================*/

typedef  cs_real_t  cs_weight_t;  /* Will allow testing single precision
                                     if set to float. */

/*============================================================================
 *  Global variables
 *============================================================================*/

bool _set[3] = {false, false, false};
cs_weight_t *_weights[3][2] = {{nullptr, nullptr},
                               {nullptr, nullptr},
                               {nullptr, nullptr}};

/* Short names for gradient computation types */

const char *cs_cell_to_vertex_type_name[]
= {N_("Unweighted"),
   N_("Shepard interpolation (weight by inverse distance)"),
   N_("Linear regression (least-squares))")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute multipliers (pseudo-weights) for unweighted case
 *
 * In this case the weights directly match the number of adjacent cells.
 *
 * \param[in]  tr_ignore  if > 0, ignore periodicity with rotation;
 *                        if > 1, ignore all periodic transforms
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_w_unweighted(int tr_ignore)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t v_sum_type = CS_DISPATCH_SUM_SIMPLE;
  if (ctx.use_gpu() == true || cs_glob_n_threads > 1)
    v_sum_type = CS_DISPATCH_SUM_ATOMIC;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0];
  int *w_sum;
  CS_REALLOC_HD(w, n_vertices, cs_weight_t, cs_alloc_mode);
  CS_MALLOC_HD(w_sum, n_vertices, int, cs_alloc_mode);

  _set[CS_CELL_TO_VERTEX_UNWEIGHTED] = true;
  _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0] = w;

  cs_array_int_fill_zero(n_vertices, w_sum);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      cs_dispatch_sum(&w_sum[v_id], 1, v_sum_type);
    }

  });

  ctx.wait();

  if (m->vtx_interfaces != nullptr)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            1,
                            true,
                            CS_INT_TYPE,
                            tr_ignore,
                            w_sum);

  ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE (cs_lnum_t v_id) {
    w[v_id] = 1. / w_sum[v_id];
  });
  ctx.wait();

  CS_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute weights based on inversed distance
 *
 * \param[in]  tr_ignore  if > 0, ignore periodicity with rotation;
 *                        if > 1, ignore all periodic transforms
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_w_inv_distance(int  tr_ignore)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t v_sum_type = CS_DISPATCH_SUM_SIMPLE;
  if (ctx.use_gpu() == true || cs_glob_n_threads > 1)
    v_sum_type = CS_DISPATCH_SUM_ATOMIC;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_real_3_t *cell_cen = mq->cell_cen;
  cs_real_t *vtx_coord = m->vtx_coord;
  const cs_real_3_t *b_face_cog = mq->b_face_cog;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_SHEPARD][0];
  cs_weight_t *wb = _weights[CS_CELL_TO_VERTEX_SHEPARD][1];
  cs_real_t *w_sum;
  CS_REALLOC_HD(w, c2v_idx[n_cells], cs_weight_t, cs_alloc_mode);
  CS_REALLOC_HD(wb, f2v_idx[n_b_faces], cs_weight_t, cs_alloc_mode);
  CS_MALLOC_HD(w_sum, n_vertices, cs_real_t, cs_alloc_mode);

  _set[CS_CELL_TO_VERTEX_SHEPARD] = true;
  _weights[CS_CELL_TO_VERTEX_SHEPARD][0] = w;
  _weights[CS_CELL_TO_VERTEX_SHEPARD][1] = wb;

  cs_array_real_fill_zero(n_vertices, w_sum);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    const cs_real_t *c_coo = cell_cen[c_id];

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      cs_real_t *v_coo = vtx_coord + v_id*3;
      cs_real_t d = cs_math_3_distance(c_coo, v_coo);
      if (d <= DBL_MIN) {
        w[j] = 1;
        w_sum[v_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1./d;
        w[j] = _w;
        cs_dispatch_sum(&w_sum[v_id], _w, v_sum_type);
      }
    }

  });

  ctx.wait();

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {

    const cs_real_t *f_coo = b_face_cog[f_id];

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = f2v_ids[j];
      cs_real_t *v_coo = vtx_coord + v_id*3;
      cs_real_t d = cs_math_3_distance(f_coo, v_coo);
      if (d <= DBL_MIN) {
        wb[j] = 1;
        w_sum[v_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1./d;
        wb[j] = _w;
        cs_dispatch_sum(&w_sum[v_id], _w, b_sum_type);
      }
    }

  });

  ctx.wait();

  if (m->vtx_interfaces != nullptr)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            1,
                            true,
                            CS_REAL_TYPE,
                            tr_ignore,
                            w_sum);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      w[j] /= w_sum[v_id];
    }

  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = f2v_ids[j];
      wb[j] /= w_sum[v_id];
    }

  }

  CS_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute factorization based on least squares
 *
 * \param[in]  tr_ignore  if > 0, ignore periodicity with rotation;
 *                        if > 1, ignore all periodic transforms
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_f_lsq(int  tr_ignore)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t v_sum_type = CS_DISPATCH_SUM_SIMPLE;
  if (ctx.use_gpu() == true || cs_glob_n_threads > 1)
    v_sum_type = CS_DISPATCH_SUM_ATOMIC;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_t *vtx_coord = m->vtx_coord;
  const cs_real_3_t *b_face_cog = mq->b_face_cog;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  cs_lnum_t  w_size = n_vertices*10;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_LR][0];
  CS_MALLOC_HD(w, w_size, cs_real_t, cs_alloc_mode);

  _set[CS_CELL_TO_VERTEX_LR] = true;
  _weights[CS_CELL_TO_VERTEX_LR][0] = w;

  cs_array_real_fill_zero(w_size, w);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t *c_coo = cell_cen[c_id];
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      const cs_real_t *v_coo = vtx_coord + v_id*3;
      cs_real_t r_coo[3]
        = {c_coo[0]-v_coo[0], c_coo[1]-v_coo[1], c_coo[2]-v_coo[2]};
      cs_real_t  *_a = w + v_id*10;
      cs_dispatch_sum(&_a[0], r_coo[0] * r_coo[0], v_sum_type); // a00
      cs_dispatch_sum(&_a[1], r_coo[1] * r_coo[0], v_sum_type); // a10
      cs_dispatch_sum(&_a[2], r_coo[1] * r_coo[1], v_sum_type); // a11
      cs_dispatch_sum(&_a[3], r_coo[2] * r_coo[0], v_sum_type); // a20
      cs_dispatch_sum(&_a[4], r_coo[2] * r_coo[1], v_sum_type); // a21
      cs_dispatch_sum(&_a[5], r_coo[2] * r_coo[2], v_sum_type); // a22
      cs_dispatch_sum(&_a[6], r_coo[0], v_sum_type);            // a30
      cs_dispatch_sum(&_a[7], r_coo[1], v_sum_type);            // a31
      cs_dispatch_sum(&_a[8], r_coo[2], v_sum_type);            // a32
      cs_dispatch_sum(&_a[9], 1., v_sum_type);                   // a33
    }
  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    const cs_real_t *f_coo = b_face_cog[f_id];
    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = f2v_ids[j];
      const cs_real_t *v_coo = vtx_coord + v_id*3;
      cs_real_t r_coo[3]
        = {f_coo[0]-v_coo[0], f_coo[1]-v_coo[1], f_coo[2]-v_coo[2]};
      cs_real_t  *_a = w + v_id*10;
      cs_dispatch_sum(&_a[0], r_coo[0] * r_coo[0], b_sum_type); // a00
      cs_dispatch_sum(&_a[1], r_coo[1] * r_coo[0], b_sum_type); // a10
      cs_dispatch_sum(&_a[2], r_coo[1] * r_coo[1], b_sum_type); // a11
      cs_dispatch_sum(&_a[3], r_coo[2] * r_coo[0], b_sum_type); // a20
      cs_dispatch_sum(&_a[4], r_coo[2] * r_coo[1], b_sum_type); // a21
      cs_dispatch_sum(&_a[5], r_coo[2] * r_coo[2], b_sum_type); // a22
      cs_dispatch_sum(&_a[6], r_coo[0], b_sum_type);            // a30
      cs_dispatch_sum(&_a[7], r_coo[1], b_sum_type);            // a31
      cs_dispatch_sum(&_a[8], r_coo[2], b_sum_type);            // a32
      cs_dispatch_sum(&_a[9], 1., b_sum_type);                   // a33
    }
  });

  ctx.wait();

  if (m->vtx_interfaces != nullptr)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            10,
                            true,
                            CS_REAL_TYPE,
                            tr_ignore,
                            w);

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    cs_math_sym_44_factor_ldlt(w + v_id*10);
}

END_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate cell values to vertex values for a strided arrray.
 *
 * \param[in]   method      interpolation method
 * \param[in]   verbosity   verbosity level
 * \param[in]   tr_ignore   if > 0, ignore periodicity with rotation;
 *                          if > 1, ignore all periodic transforms
 * \param[in]   c_weight    cell weight, or nullptr
 * \param[in]   c_var       base cell-based variable
 * \param[in]   b_var       base boundary-face values, or nullptr
 * \param[out]  v_var       vertex-based variable
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_cell_to_vertex_strided(cs_cell_to_vertex_type_t   method,
                        int                        verbosity,
                        int                        tr_ignore,
                        const cs_real_t *restrict  c_weight,
                        const cs_real_t *restrict  c_var,
                        const cs_real_t *restrict  b_var,
                        cs_real_t *restrict        v_var)
{
  CS_UNUSED(verbosity);

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t v_sum_type = CS_DISPATCH_SUM_SIMPLE;
  if (ctx.use_gpu() == true || cs_glob_n_threads > 1)
    v_sum_type = CS_DISPATCH_SUM_ATOMIC;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict vtx_coord = m->vtx_coord;
  const cs_real_3_t *restrict b_face_cog = mq->b_face_cog;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  const cs_lnum_t n_v_values = n_vertices * stride;

  const cs_real_3_t *cell_cen = mq->cell_cen;

  cs_array_real_fill_zero(n_v_values, v_var);

  switch(method) {

  case CS_CELL_TO_VERTEX_UNWEIGHTED:
    {
      if (c_weight == nullptr) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id*stride + k], c_var[c_id*stride + k], v_sum_type);
          }
        });
        ctx.wait();

        if (m->vtx_interfaces != nullptr)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  m->n_vertices,
                                  stride,
                                  true,
                                  CS_REAL_TYPE,
                                  tr_ignore,
                                  v_var);

        if (! _set[CS_CELL_TO_VERTEX_UNWEIGHTED])
          _cell_to_vertex_w_unweighted(tr_ignore);

        const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0];
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE (cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id*stride + k] *= w[v_id];
        });
      }
      else {
        cs_real_t *v_w;
        CS_MALLOC(v_w, n_vertices, cs_real_t);
        cs_array_real_fill_zero(n_vertices, v_w);
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id*stride + k],
                              c_var[c_id*stride + k] * c_weight[c_id],
                              v_sum_type);
            cs_dispatch_sum(&v_w[v_id], c_weight[c_id], v_sum_type);
          }
        });
        ctx.wait();

        if (m->vtx_interfaces != nullptr) {
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  stride,
                                  true,
                                  CS_REAL_TYPE,
                                  tr_ignore,
                                  v_var);
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  tr_ignore,
                                  v_w);
        }
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE (cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id*stride + k] /= v_w[v_id];
        });
        ctx.wait();
        CS_FREE(v_w);
      }
    }
    break;

  case CS_CELL_TO_VERTEX_SHEPARD:
    {
      if (! _set[CS_CELL_TO_VERTEX_SHEPARD])
        _cell_to_vertex_w_inv_distance(tr_ignore);

      const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_SHEPARD][0];

      cs_real_t *v_w = nullptr;
      if (c_weight != nullptr) {
        CS_MALLOC(v_w, n_vertices, cs_real_t);
        cs_array_real_fill_zero(n_vertices, v_w);
      }

      if (c_weight == nullptr) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id*stride + k],
                              c_var[c_id*stride + k] * w[j],
                              v_sum_type);
          }
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id*stride + k],
                              c_var[c_id*stride + k] * w[j] * c_weight[c_id],
                              v_sum_type);
            cs_dispatch_sum(&v_w[v_id], w[j] * c_weight[c_id], v_sum_type);
          }
        });
      }

      const cs_weight_t *wb = _weights[CS_CELL_TO_VERTEX_SHEPARD][1];

      if (c_weight == nullptr) {

        if (b_var == nullptr) {
          ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
            cs_lnum_t c_id = b_face_cells[f_id];
            const cs_real_t *_b_var = c_var + c_id*stride;
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < stride; k++)
                cs_dispatch_sum(&v_var[v_id*stride+k], _b_var[k] * wb[j], b_sum_type);
            }
          });
        }
        else {
          ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < stride; k++)
                cs_dispatch_sum(&v_var[v_id*stride+k],
                                b_var[f_id*stride+k] * wb[j],
                                b_sum_type);
            }
          });
        }
      }
      else { /* c_weight != nullptr */

        if (b_var == nullptr) {
          ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            const cs_real_t *_b_var = c_var + c_id*stride;
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < stride; k++)
                cs_dispatch_sum(&v_var[v_id*stride+k],
                                _b_var[k] * wb[j] * c_weight[c_id],
                                b_sum_type);
              cs_dispatch_sum(&v_w[v_id], wb[j] * c_weight[c_id], b_sum_type);
            }
          });
        }
        else {
          ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < stride; k++)
                cs_dispatch_sum(&v_var[v_id*stride+k],
                                b_var[f_id*stride+k] * wb[j] * c_weight[c_id],
                                b_sum_type);
              cs_dispatch_sum(&v_w[v_id], wb[j] * c_weight[c_id], b_sum_type);
            }
          });
        }
      }

      ctx.wait();

      if (m->vtx_interfaces != nullptr)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                stride,
                                true,
                                CS_REAL_TYPE,
                                tr_ignore,
                                v_var);

      if (c_weight != nullptr) {
        if (m->vtx_interfaces != nullptr)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  tr_ignore,
                                  v_w);
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE (cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id*stride+k] /= v_w[v_id];
        });

        ctx.wait();
        CS_FREE(v_w);
      }

    }
    break;

  case CS_CELL_TO_VERTEX_LR:
    {
      if (! _set[CS_CELL_TO_VERTEX_LR])
        _cell_to_vertex_f_lsq(tr_ignore);

      cs_real_t  *rhs;
      cs_lnum_t  rhs_size = n_vertices*4*stride;
      CS_MALLOC_HD(rhs, rhs_size, cs_real_t, cs_alloc_mode);
      cs_array_real_fill_zero(rhs_size, rhs);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        const cs_real_t *c_coo = cell_cen[c_id];
        const cs_real_t *_c_var = c_var + c_id*stride;
        cs_lnum_t s_id = c2v_idx[c_id];
        cs_lnum_t e_id = c2v_idx[c_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t v_id = c2v_ids[j];
          const cs_real_t *v_coo = vtx_coord + v_id*3;
          cs_real_t r_coo[3]
            = {c_coo[0]-v_coo[0], c_coo[1]-v_coo[1], c_coo[2]-v_coo[2]};
          for (cs_lnum_t k = 0; k < stride; k++) {
            cs_real_t *_rhs = rhs + v_id*4*stride + k*4;
            cs_dispatch_sum(&_rhs[0], r_coo[0] * _c_var[k], v_sum_type);
            cs_dispatch_sum(&_rhs[1], r_coo[1] * _c_var[k], v_sum_type);
            cs_dispatch_sum(&_rhs[2], r_coo[2] * _c_var[k], v_sum_type);
            cs_dispatch_sum(&_rhs[3], _c_var[k], v_sum_type);
          }
        }
      });

      if (b_var == nullptr) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
          const cs_real_t *f_coo = b_face_cog[f_id];
          cs_lnum_t c_id = b_face_cells[f_id];
          const cs_real_t *_b_var = c_var + c_id*stride;
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            const cs_real_t *v_coo = vtx_coord + v_id*3;
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_real_t *_rhs = rhs + v_id*4*stride + k*4;
              cs_dispatch_sum(&_rhs[0], (f_coo[0] - v_coo[0]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[1], (f_coo[1] - v_coo[1]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[2], (f_coo[2] - v_coo[2]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[3], _b_var[k], b_sum_type);
            }
          }
        });
      }
      else {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
          const cs_real_t *f_coo = b_face_cog[f_id];
          const cs_real_t *_b_var = b_var + f_id*stride;
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            const cs_real_t *v_coo = vtx_coord + v_id*3;
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_real_t *_rhs = rhs + v_id*4*stride + k*4;
              cs_dispatch_sum(&_rhs[0], (f_coo[0] - v_coo[0]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[1], (f_coo[1] - v_coo[1]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[2], (f_coo[2] - v_coo[2]) * _b_var[k], b_sum_type);
              cs_dispatch_sum(&_rhs[3], _b_var[k], b_sum_type);
            }
          }
        });
      }

      ctx.wait();

      if (m->vtx_interfaces != nullptr)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                4*stride,
                                true,
                                CS_REAL_TYPE,
                                tr_ignore,
                                rhs);

      const cs_weight_t *ldlt = _weights[CS_CELL_TO_VERTEX_LR][0];

      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
        const cs_real_t *_ldlt = ldlt + v_id*10;
        for (cs_lnum_t k = 0; k < stride; k++) {
          const cs_real_t *_rhs = rhs + v_id*4*stride + k*4;
          v_var[v_id*stride + k]
            = cs_math_sym_44_partial_solve_ldlt(_ldlt, _rhs);
        }
      }

      ctx.wait();
      CS_FREE(rhs);
    }
    break;

  default:
    break;
  }

  ctx.wait();
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free cell to vertex interpolation weights.
 *
 * This will force subsequent calls to rebuild those weights if needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_to_vertex_free(void)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++)
    CS_FREE(_weights[i][j]);
    _set[i] = false;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate cell values to vertex values.
 *
 * \param[in]       method            interpolation method
 * \param[in]       verbosity         verbosity level
 * \param[in]       ignore_rot_perio  if true, ignore periodicity of rotation
 * \param[in]       c_weight          cell weight, or nullptr
 * \param[in]       c_var             base cell-based variable
 * \param[in]       b_var             base boundary-face values, or nullptr
 * \param[out]      v_var             vertex-based variable
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_cell_to_vertex(cs_cell_to_vertex_type_t   method,
                  int                        verbosity,
                  bool                       ignore_rot_perio,
                  const cs_real_t *restrict  c_weight,
                  const cs_real_t *restrict  c_var,
                  const cs_real_t *restrict  b_var,
                  cs_real_t *restrict        v_var)
{
  CS_UNUSED(verbosity);

  int tr_ignore = (ignore_rot_perio) ? 1 : 0;

  _cell_to_vertex_strided<stride>(method,
                                  verbosity,
                                  tr_ignore,
                                  c_weight,
                                  c_var,
                                  b_var,
                                  v_var);
}

// Force instanciation

template void
cs_cell_to_vertex<1>(cs_cell_to_vertex_type_t  method,
                  int                          verbosity,
                  bool                         ignore_rot_perio,
                  const cs_real_t *restrict    c_weight,
                  const cs_real_t *restrict    c_var,
                  const cs_real_t *restrict    b_var,
                  cs_real_t *restrict          v_var);

template void
cs_cell_to_vertex<3>(cs_cell_to_vertex_type_t  method,
                  int                          verbosity,
                  bool                         ignore_rot_perio,
                  const cs_real_t *restrict    c_weight,
                  const cs_real_t *restrict    c_var,
                  const cs_real_t *restrict    b_var,
                  cs_real_t *restrict          v_var);

template void
cs_cell_to_vertex<6>(cs_cell_to_vertex_type_t  method,
                  int                          verbosity,
                  bool                         ignore_rot_perio,
                  const cs_real_t *restrict    c_weight,
                  const cs_real_t *restrict    c_var,
                  const cs_real_t *restrict    b_var,
                  cs_real_t *restrict          v_var);

/*----------------------------------------------------------------------------*/
