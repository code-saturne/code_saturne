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
 * Standard library headers
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

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for gradient computation types */

const char *cs_face_to_vertex_type_name[] = {
  N_("Unweighted"),
  N_("Surface barycenter (weight by surface)"),
  N_("Shepard interpolation (weight by inverse distance)"),
};

/*============================================================================
 * Public function definitions
 *============================================================================*/

template <cs_lnum_t stride>
void
cs_face_to_vertex_t<stride>::initialize(cs_lnum_t        n_faces,
                                        const cs_lnum_t *list_faces,
                                        cs_lnum_t        n_vtx,
                                        const cs_lnum_t *list_vtx)
{
  _n_faces    = n_faces;
  _list_faces = list_faces;
  _n_vtx      = n_vtx;
  _list_vtx   = list_vtx;

  const cs_mesh_t *m          = cs_glob_mesh;
  const cs_lnum_t  n_vertices = m->n_vertices;

  if (_list_vtx == nullptr && n_vertices != _n_vtx) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: %s incompatible number of verticies."),
              __func__);
  }

  if (_list_faces == nullptr && m->n_b_faces != _n_faces) {
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: %s incompatible number of boundary faces."),
              __func__);
  }

  CS_MALLOC(_v_w, n_vertices, cs_real_t);
  CS_MALLOC(_v_var, stride * n_vertices, cs_real_t);

  cs_dispatch_context ctx;

  ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
    _v_w[v_id] = 0.;
    for (cs_lnum_t k = 0; k < stride; ++k) {
      _v_var[v_id * stride + k] = 0.0;
    }
  });
  ctx.wait();
};

template <cs_lnum_t stride>
void
cs_face_to_vertex_t<stride>::free()
{
  CS_FREE(_v_w);
  CS_FREE(_v_var);
};

template <cs_lnum_t stride>
void
cs_face_to_vertex_t<stride>::compute_on_boundary(
  cs_face_to_vertex_type_t method,
  bool                     ignore_rot_perio,
  const cs_real_t         *b_var,
  cs_real_t               *v_var)
{
  const cs_mesh_t            *m        = cs_glob_mesh;
  const cs_mesh_quantities_t *mq       = cs_glob_mesh_quantities;
  const cs_real_t            *vtx_coor = m->vtx_coord;

  cs_dispatch_context    ctx;
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_lnum_t n_vertices = m->n_vertices;

  const cs_lnum_t *bf2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *bf2v_ids = m->b_face_vtx_lst;

  // Initialize
  ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
    _v_w[v_id] = 0.;
    for (cs_lnum_t k = 0; k < stride; ++k) {
      _v_var[v_id * stride + k] = 0.0;
    }
  });
  ctx.wait();

  switch (method) {
    case CS_FACE_TO_VERTEX_UNWEIGHTED: {
      if (_list_faces == nullptr && _n_faces > 0) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_id * stride + k],
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], 1.0, b_sum_type);
          }
        });
      }
      else {
        ctx.parallel_for(_n_faces, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_l_id) {
          cs_lnum_t bf_id = _list_faces[bf_l_id];
          cs_lnum_t s_id  = bf2v_idx[bf_id];
          cs_lnum_t e_id  = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_l_id * stride + k],
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], 1.0, b_sum_type);
          }
        });
      }

      ctx.wait();
    } break;

    case CS_FACE_TO_VERTEX_SURFACE: {
      const cs_real_t *bf_surf = mq->b_face_surf;

      if (_list_faces == nullptr && _n_faces > 0) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          cs_real_t bw   = bf_surf[bf_id];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * bw,
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], bw, b_sum_type);
          }
        });
      }
      else {
        ctx.parallel_for(_n_faces, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_l_id) {
          cs_lnum_t bf_id = _list_faces[bf_l_id];
          cs_lnum_t s_id  = bf2v_idx[bf_id];
          cs_lnum_t e_id  = bf2v_idx[bf_id + 1];
          cs_real_t bw    = bf_surf[bf_id];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = bf2v_ids[j];
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_l_id * stride + k] * bw,
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], bw, b_sum_type);
          }
        });
      }

      ctx.wait();
    } break;

    case CS_FACE_TO_VERTEX_SHEPARD: {
      const cs_real_3_t *b_face_cog = mq->b_face_cog;

      if (_list_faces == nullptr && _n_faces > 0) {
        ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_id) {
          const cs_real_t *bf_coor = b_face_cog[bf_id];

          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t        v_id   = bf2v_ids[j];
            const cs_real_t *v_coor = vtx_coor + v_id * 3;
            cs_real_t        dist   = cs_math_3_distance(bf_coor, v_coor);
            // Power 1 as cell_to_vertex
            cs_real_t bw = 1.0 / dist;
            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_id * stride + k] * bw,
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], bw, b_sum_type);
          }
        });
      }
      else {
        ctx.parallel_for(_n_faces, [=] CS_F_HOST_DEVICE(cs_lnum_t bf_l_id) {
          cs_lnum_t        bf_id   = _list_faces[bf_l_id];
          const cs_real_t *bf_coor = b_face_cog[bf_id];

          cs_lnum_t s_id = bf2v_idx[bf_id];
          cs_lnum_t e_id = bf2v_idx[bf_id + 1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t        v_id   = bf2v_ids[j];
            const cs_real_t *v_coor = vtx_coor + v_id * 3;
            cs_real_t        dist   = cs_math_3_distance(bf_coor, v_coor);
            // Power 1 as cell_to_vertex
            cs_real_t bw = 1.0 / dist;

            for (cs_lnum_t k = 0; k < stride; k++) {
              cs_dispatch_sum(&_v_var[v_id * stride + k],
                              b_var[bf_l_id * stride + k] * bw,
                              b_sum_type);
            }
            cs_dispatch_sum(&_v_w[v_id], bw, b_sum_type);
          }
        });
      }

      ctx.wait();
    } break;

    default:
      bft_error(__FILE__,
                __LINE__,
                0,
                _("Error: %s option not available."),
                __func__);
      break;
  }

  if (m->vtx_interfaces != nullptr) {
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            n_vertices,
                            stride,
                            true,
                            CS_REAL_TYPE,
                            ignore_rot_perio,
                            _v_var);

    cs_interface_set_sum_tr(m->vtx_interfaces,
                            n_vertices,
                            1,
                            true,
                            CS_REAL_TYPE,
                            ignore_rot_perio,
                            _v_w);
  }

  if (_list_vtx == nullptr && _n_vtx > 0) {
    ctx.parallel_for(n_vertices, [=] CS_F_HOST_DEVICE(cs_lnum_t v_id) {
      cs_real_t ivw = 1.0;
      if (_v_w[v_id] > cs_dbl_epsilon) {
        ivw = 1.0 / _v_w[v_id];
      }
      for (cs_lnum_t k = 0; k < stride; k++) {
        v_var[v_id * stride + k] = _v_var[v_id * stride + k] * ivw;
      }
    });
  }
  else {
    ctx.parallel_for(_n_vtx, [=] CS_F_HOST_DEVICE(cs_lnum_t v_l_id) {
      cs_lnum_t v_id = _list_vtx[v_l_id];
      cs_real_t ivw  = 1.0;
      if (_v_w[v_id] > cs_dbl_epsilon) {
        ivw = 1.0 / _v_w[v_id];
      }
      for (cs_lnum_t k = 0; k < stride; k++) {
        v_var[v_l_id * stride + k] = _v_var[v_id * stride + k] * ivw;
      }
    });
  }

  ctx.wait();
};

template class cs_face_to_vertex_t<3>;

/*----------------------------------------------------------------------------*/
