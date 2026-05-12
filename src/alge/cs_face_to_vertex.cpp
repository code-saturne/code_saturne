/*============================================================================
 * Face to vertex interpolation.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

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
#include "base/cs_timer.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_face_to_vertex.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_face_to_vertex.cpp
 * \brief Cell to vertex interpolation.
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL (CS_CL_SIZE / 8)

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.E+12
#endif

/*=============================================================================
 * Local type definitions
 *============================================================================*/

typedef cs_real_t cs_weight_t; /* Will allow testing single precision
                                  if set to float. */

/*============================================================================
 *  Global variables
 *============================================================================*/

bool         _set_bface[3]     = { false, false, false };
cs_weight_t *_weights_bface[3] = { nullptr, nullptr, nullptr };

/* Short names for gradient computation types */

const char *cs_face_to_vertex_type_name[] = {
  N_("Unweighted"),
  N_("Surface barycenter (weight by surface)"),
  N_("Shepard interpolation (weight by inverse distance)"),
};

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
_bface_to_vertex_w_unweighted(int tr_ignore)
{
  const cs_mesh_t *m = cs_glob_mesh;

  cs_dispatch_context    ctx;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_lnum_t n_vertices = m->n_vertices;

  const cs_lnum_t *bf2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *bf2v_ids = m->b_face_vtx_lst;

  cs_weight_t *w = _weights_bface[CS_FACE_TO_VERTEX_UNWEIGHTED];
  int         *w_sum;
  CS_REALLOC_HD(w, n_vertices, cs_weight_t, cs_alloc_mode);
  CS_MALLOC_HD(w_sum, n_vertices, int, cs_alloc_mode);

  _set_bface[CS_FACE_TO_VERTEX_UNWEIGHTED]     = true;
  _weights_bface[CS_FACE_TO_VERTEX_UNWEIGHTED] = w;

  cs_array_int_fill_zero(n_vertices, w_sum);

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
    cs_lnum_t s_id = bf2v_idx[bf_id];
    cs_lnum_t e_id = bf2v_idx[bf_id + 1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = bf2v_ids[j];
      cs_dispatch_sum(&w_sum[v_id], 1, b_sum_type);
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

  ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
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
_bface_to_vertex_w_inv_distance(int tr_ignore)
{
  const cs_mesh_t            *m  = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_dispatch_context    ctx;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  cs_real_t         *vtx_coord  = m->vtx_coord;
  const cs_real_3_t *b_face_cog = mq->b_face_cog;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_b_faces  = m->n_b_faces;

  const cs_lnum_t *bf2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *bf2v_ids = m->b_face_vtx_lst;

  cs_weight_t *wb = _weights_bface[CS_FACE_TO_VERTEX_SHEPARD];
  cs_real_t   *w_sum;
  CS_REALLOC_HD(wb, bf2v_idx[n_b_faces], cs_weight_t, cs_alloc_mode);
  CS_MALLOC_HD(w_sum, n_vertices, cs_real_t, cs_alloc_mode);

  _set_bface[CS_FACE_TO_VERTEX_SHEPARD]     = true;
  _weights_bface[CS_FACE_TO_VERTEX_SHEPARD] = wb;

  ctx.parallel_for(n_vertices,
                   [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) { w_sum[v_id] = 0.; });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
    const cs_real_t *f_coo = b_face_cog[bf_id];

    cs_lnum_t s_id = bf2v_idx[bf_id];
    cs_lnum_t e_id = bf2v_idx[bf_id + 1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t  v_id  = bf2v_ids[j];
      cs_real_t *v_coo = vtx_coord + v_id * 3;
      cs_real_t  d     = cs_math_3_distance(f_coo, v_coo);
      if (d <= DBL_MIN) {
        wb[j]       = 1;
        w_sum[v_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1. / d;
        wb[j]        = _w;
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

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
    cs_lnum_t s_id = bf2v_idx[bf_id];
    cs_lnum_t e_id = bf2v_idx[bf_id + 1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = bf2v_ids[j];
      wb[j] /= w_sum[v_id];
    }
  });

  ctx.wait();

  CS_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate boundary face values to vertex values for a strided
 * arrray.
 *
 * \param[in]   method      interpolation method
 * \param[in]   verbosity   verbosity level
 * \param[in]   tr_ignore   if > 0, ignore periodicity with rotation;
 *                          if > 1, ignore all periodic transforms
 * \param[in]   b_weight    boundary-face weight, or nullptr
 * \param[in]   b_var       base boundary-face values
 * \param[out]  v_var       vertex-based variable
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_bface_to_vertex_strided(cs_face_to_vertex_type_t method,
                         [[maybe_unused]] int     verbosity,
                         int                      tr_ignore,
                         const cs_real_t *restrict b_weight,
                         const cs_real_t *restrict b_var,
                         cs_real_t *restrict v_var)
{
  const cs_mesh_t            *m  = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_dispatch_context    ctx;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_lnum_t n_vertices = m->n_vertices;

  const cs_lnum_t *bf2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *bf2v_ids = m->b_face_vtx_lst;

  const cs_lnum_t n_v_values = n_vertices * stride;

  ctx.parallel_for(n_v_values,
                   [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) { v_var[v_id] = 0.; });
  ctx.wait();

  switch (method) {
    case CS_FACE_TO_VERTEX_UNWEIGHTED: {
      if (b_weight == nullptr) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k],
                              b_sum_type);
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

        if (!_set_bface[CS_FACE_TO_VERTEX_UNWEIGHTED])
          _bface_to_vertex_w_unweighted(tr_ignore);

        const cs_weight_t *wb = _weights_bface[CS_FACE_TO_VERTEX_UNWEIGHTED];
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id * stride + k] *= wb[v_id];
        });
        ctx.wait();
      }
      else {
        cs_real_t *v_w;
        CS_MALLOC(v_w, n_vertices, cs_real_t);

        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
          v_w[v_id] = 0.;
        });

        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * b_weight[bf_id],
                              b_sum_type);
            cs_dispatch_sum(&v_w[v_id], b_weight[bf_id], b_sum_type);
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
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id * stride + k] /= v_w[v_id];
        });
        ctx.wait();
        CS_FREE(v_w);
      }
    } break;

    case CS_FACE_TO_VERTEX_SURFACE: {
      const cs_real_t *bf_surf = mq->b_face_surf;

      cs_real_t *v_w = nullptr;
      CS_MALLOC(v_w, n_vertices, cs_real_t);

      ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
        v_w[v_id] = 0.;
      });
      ctx.wait();

      if (b_weight == nullptr) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * bf_surf[bf_id],
                              b_sum_type);
            }
            cs_dispatch_sum(&v_w[v_id], bf_surf[bf_id], b_sum_type);
          }
        });
      }
      else { /* b_weight != nullptr */

        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * bf_surf[bf_id] *
                                b_weight[bf_id],
                              b_sum_type);
            }
            cs_dispatch_sum(&v_w[v_id],
                            bf_surf[bf_id] * b_weight[bf_id],
                            b_sum_type);
          }
        });
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

      if (m->vtx_interfaces != nullptr)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                n_vertices,
                                1,
                                true,
                                CS_REAL_TYPE,
                                tr_ignore,
                                v_w);
      ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
        for (cs_lnum_t k = 0; k < stride; k++)
          v_var[v_id * stride + k] /= v_w[v_id];
      });

      ctx.wait();
      CS_FREE(v_w);

    } break;

    case CS_FACE_TO_VERTEX_SHEPARD: {
      if (!_set_bface[CS_FACE_TO_VERTEX_SHEPARD])
        _bface_to_vertex_w_inv_distance(tr_ignore);

      const cs_weight_t *wb = _weights_bface[CS_FACE_TO_VERTEX_SHEPARD];

      cs_real_t *v_w = nullptr;
      if (b_weight != nullptr) {
        CS_MALLOC(v_w, n_vertices, cs_real_t);

        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
          v_w[v_id] = 0.;
        });
      }

      if (b_weight == nullptr) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * wb[j],
                              b_sum_type);
          }
        });
      }
      else { /* b_weight != nullptr */

        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++)
              cs_dispatch_sum(&v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * wb[j] *
                                b_weight[bf_id],
                              b_sum_type);
            cs_dispatch_sum(&v_w[v_id], wb[j] * b_weight[bf_id], b_sum_type);
          }
        });
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

      if (b_weight != nullptr) {
        if (m->vtx_interfaces != nullptr)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  tr_ignore,
                                  v_w);
        ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
          for (cs_lnum_t k = 0; k < stride; k++)
            v_var[v_id * stride + k] /= v_w[v_id];
        });

        ctx.wait();
        CS_FREE(v_w);
      }

    } break;
    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                _("Error: %s option not available."),
                __func__);
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
 * \brief  Free face to vertex interpolation weights.
 *
 * This will force subsequent calls to rebuild those weights if needed.
 */
/*----------------------------------------------------------------------------*/

void
cs_face_to_vertex_free(void)
{
  for (int i = 0; i < 3; i++) {
    CS_FREE(_weights_bface[i]);
    _set_bface[i] = false;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate boundary face values to vertex values.
 *
 * \param[in]       method            interpolation method
 * \param[in]       verbosity         verbosity level
 * \param[in]       ignore_rot_perio  if true, ignore periodicity of rotation
 * \param[in]       b_weight          boundary-face weight, or nullptr
 * \param[in]       b_var             base boundary-face values, or nullptr
 * \param[out]      v_var             vertex-based variable
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_bface_to_vertex(cs_face_to_vertex_type_t method,
                   int                      verbosity,
                   bool                     ignore_rot_perio,
                   const cs_real_t *restrict b_weight,
                   const cs_real_t *restrict b_var,
                   cs_real_t *restrict v_var)
{
  int tr_ignore = (ignore_rot_perio) ? 1 : 0;

  _bface_to_vertex_strided<stride>(method,
                                   verbosity,
                                   tr_ignore,
                                   b_weight,
                                   b_var,
                                   v_var);
}

// Force instanciation

template void
cs_bface_to_vertex<1>(cs_face_to_vertex_type_t method,
                      int                      verbosity,
                      bool                     ignore_rot_perio,
                      const cs_real_t *restrict b_weight,
                      const cs_real_t *restrict b_var,
                      cs_real_t *restrict v_var);

template void
cs_bface_to_vertex<3>(cs_face_to_vertex_type_t method,
                      int                      verbosity,
                      bool                     ignore_rot_perio,
                      const cs_real_t *restrict b_weight,
                      const cs_real_t *restrict b_var,
                      cs_real_t *restrict v_var);

template void
cs_bface_to_vertex<6>(cs_face_to_vertex_type_t method,
                      int                      verbosity,
                      bool                     ignore_rot_perio,
                      const cs_real_t *restrict b_weight,
                      const cs_real_t *restrict b_var,
                      cs_real_t *restrict v_var);

/*----------------------------------------------------------------------------*/
