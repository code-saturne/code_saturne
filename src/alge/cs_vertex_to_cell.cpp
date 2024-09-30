/*============================================================================
 * Vertex to cell interpolation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_vertex_to_cell.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_vertex_to_cell.c
 * \brief Cell to vertex interpolation..
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/*=============================================================================
 * Local type definitions
 *============================================================================*/

typedef  cs_real_t  cs_weight_t;  /* will allow testing single precision
                                     if set to float */

/*============================================================================
 *  Global variables
 *============================================================================*/

bool          _set_vtc[3] = {false, false, false};
cs_weight_t  *_weights_vtc[3] = {nullptr, nullptr, nullptr};

/* Short names for gradient computation types */

const char *cs_vertex_to_cell_type_name[]
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
 * In this case the weights directly match the number of adjacent vertices.
 */
/*----------------------------------------------------------------------------*/

static void
_vertex_to_cell_w_unweighted(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;

  cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED];
  BFT_REALLOC(w, n_cells, cs_weight_t);

  _set_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED] = true;
  _weights_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED] = w;

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    w[c_id] = 1. / (c2v_idx[c_id+1]-c2v_idx[c_id]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute weights based on inversed distance
 */
/*----------------------------------------------------------------------------*/

static void
_vertex_to_cell_w_inv_distance(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_SHEPARD];
  cs_real_t   *w_sum;
  BFT_REALLOC(w, c2v_idx[n_cells], cs_weight_t);
  BFT_MALLOC(w_sum, n_cells, cs_real_t);

  _set_vtc[CS_VERTEX_TO_CELL_SHEPARD] = true;
  _weights_vtc[CS_VERTEX_TO_CELL_SHEPARD] = w;

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    w_sum[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t *c_coo = mq->cell_cen + c_id*3;

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      cs_real_t *v_coo = m->vtx_coord + v_id*3;
      cs_real_t d = cs_math_3_distance(v_coo, c_coo);
      if (d <= DBL_MIN) {
        w[j] = 1;
        w_sum[c_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1./d;
        w[j] = _w;
        w_sum[c_id] += _w;
      }
    }

  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      w[j] /= w_sum[c_id];
    }

  }

  BFT_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute factorization based on least squares
 */
/*----------------------------------------------------------------------------*/

static void
_vertex_to_cell_f_lsq(void)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  cs_lnum_t  w_size = n_cells*10;

  cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_LR];
  BFT_MALLOC(w, w_size, cs_real_t);

  _set_vtc[CS_VERTEX_TO_CELL_LR] = true;
  _weights_vtc[CS_VERTEX_TO_CELL_LR] = w;

# pragma omp parallel for if(w_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < w_size; i++)
    w[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t *c_coo = mq->cell_cen + c_id*3;
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      const cs_real_t *v_coo = m->vtx_coord + v_id*3;
      cs_real_t r_coo[3]
        = {v_coo[0]-c_coo[0], v_coo[1]-c_coo[1], v_coo[2]-c_coo[2]};
      cs_real_t  *_a = w + c_id*10;
      _a[0] += r_coo[0] * r_coo[0]; // a00
      _a[1] += r_coo[1] * r_coo[0]; // a10
      _a[2] += r_coo[1] * r_coo[1]; // a11
      _a[3] += r_coo[2] * r_coo[0]; // a20
      _a[4] += r_coo[2] * r_coo[1]; // a21
      _a[5] += r_coo[2] * r_coo[2]; // a22
      _a[6] += r_coo[0];            // a30
      _a[7] += r_coo[1];            // a31
      _a[8] += r_coo[2];            // a32
      _a[9] += 1;                   // a33
    }
  }

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    cs_math_sym_44_factor_ldlt(w + c_id*10);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate vertex values to cell values for a scalar arrray.
 *
 * \param[in]  method       interpolation method
 * \param[in]  verbosity    verbosity level
 * \param[in]  v_weight     vertex weight, or nullptr
 * \param[in]  v_var        base vertex-based variable
 * \param[out] c_var        cell-based variable
 */
/*----------------------------------------------------------------------------*/

static void
_vertex_to_cell_scalar(cs_vertex_to_cell_type_t   method,
                       int                        verbosity,
                       const cs_real_t *restrict  v_weight,
                       const cs_real_t *restrict  v_var,
                       cs_real_t *restrict        c_var)
{
  CS_UNUSED(verbosity);

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_mesh_adjacencies_cell_vertices();

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    c_var[c_id] = 0;

  switch(method) {

  case CS_VERTEX_TO_CELL_UNWEIGHTED:
    {
      if (v_weight == nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            c_var[c_id] += v_var[v_id];
          }
        }

        if (! _set_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED])
          _vertex_to_cell_w_unweighted();

        const cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED];
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_var[c_id] *= w[c_id];
      }
      else {
        cs_real_t *c_w;
        BFT_MALLOC(c_w, n_cells, cs_real_t);
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          c_w[c_id] = 0;
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            c_var[c_id] += v_var[v_id]*v_weight[v_id];
            c_w[c_id] += v_weight[v_id];
          }
        }

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_var[c_id] /= c_w[c_id];

        BFT_FREE(c_w);
      }
    }
    break;

  case CS_VERTEX_TO_CELL_SHEPARD:
    {
      if (! _set_vtc[CS_VERTEX_TO_CELL_SHEPARD])
        _vertex_to_cell_w_inv_distance();

      const cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_SHEPARD];

      cs_real_t *c_w = nullptr;
      if (v_weight != nullptr) {
        BFT_MALLOC(c_w, n_cells, cs_real_t);
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_w[c_id] = 0;
      }

      if (v_weight == nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            c_var[c_id] += v_var[v_id] * w[j];
          }
        }
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            c_var[c_id] += v_var[v_id] * w[j] * v_weight[v_id];
            c_w[c_id] += w[j] * v_weight[v_id];
          }
        }
      }

      if (v_weight != nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_var[c_id] /= c_w[c_id];
        BFT_FREE(c_w);
      }

    }
    break;

  case CS_VERTEX_TO_CELL_LR:
    {
      if (! _set_vtc[CS_VERTEX_TO_CELL_LR])
        _vertex_to_cell_f_lsq();

      cs_real_t  *rhs;
      cs_lnum_t  rhs_size = n_cells*4;
      BFT_MALLOC(rhs, rhs_size, cs_real_t);
      for (cs_lnum_t i = 0; i < rhs_size; i++)
        rhs[i] = 0;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t *c_coo = mq->cell_cen + c_id*3;
        cs_lnum_t s_id = c2v_idx[c_id];
        cs_lnum_t e_id = c2v_idx[c_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t v_id = c2v_ids[j];
          cs_real_t _v_var = v_var[v_id];
          const cs_real_t *v_coo = m->vtx_coord + v_id*3;
          cs_real_t r_coo[3]
            = {v_coo[0]-c_coo[0], v_coo[1]-c_coo[1], v_coo[2]-c_coo[2]};
          cs_real_t  *_rhs = rhs + c_id*4;
          _rhs[0] += r_coo[0] * _v_var;
          _rhs[1] += r_coo[1] * _v_var;
          _rhs[2] += r_coo[2] * _v_var;
          _rhs[3] += _v_var;
        }
      }

      const cs_weight_t *ldlt = _weights_vtc[CS_VERTEX_TO_CELL_LR];

#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t  *_ldlt = ldlt + c_id*10;
        const cs_real_t  *_rhs = rhs + c_id*4;
        c_var[c_id] = cs_math_sym_44_partial_solve_ldlt(_ldlt, _rhs);
      }

      BFT_FREE(rhs);
    }
    break;
  default:
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate vertex values to cell values for a strided arrray.
 *
 * \param[in]   method      interpolation method
 * \param[in]   verbosity   verbosity level
 * \param[in]   var_dim     variable dimension
 * \param[in]   v_weight    vertex weight, or nullptr
 * \param[in]   v_var       base vertex-based variable
 * \param[out]  c_var       cell-based variable
 */
/*----------------------------------------------------------------------------*/

static void
_vertex_to_cell_strided(cs_vertex_to_cell_type_t   method,
                        int                        verbosity,
                        cs_lnum_t                  var_dim,
                        const cs_real_t *restrict  v_weight,
                        const cs_real_t *restrict  v_var,
                        cs_real_t *restrict        c_var)
{
  CS_UNUSED(verbosity);

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_mesh_adjacencies_cell_vertices();

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t n_c_values = n_cells*var_dim;

# pragma omp parallel for if(n_c_values > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_c_values; c_id++)
    c_var[c_id] = 0;

  switch(method) {

  case CS_VERTEX_TO_CELL_UNWEIGHTED:
    {
      if (v_weight == nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              c_var[c_id*var_dim + k] += v_var[v_id*var_dim + k];
          }
        }

        if (! _set_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED])
          _vertex_to_cell_w_unweighted();

        const cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_UNWEIGHTED];
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t k = 0; k < var_dim; k++)
            c_var[c_id*var_dim + k] *= w[c_id];
        }
      }
      else {
        cs_real_t *c_w;
        BFT_MALLOC(c_w, n_cells, cs_real_t);
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_w[c_id] = 0;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              c_var[c_id*var_dim + k] += v_var[v_id*var_dim + k] * v_weight[v_id];
            c_w[c_id] += v_weight[v_id];
          }
        }

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t k = 0; k < var_dim; k++)
            c_var[c_id*var_dim + k] /= c_w[c_id];
        }

        BFT_FREE(c_w);
      }
    }
    break;

  case CS_VERTEX_TO_CELL_SHEPARD:
    {
      if (! _set_vtc[CS_VERTEX_TO_CELL_SHEPARD])
        _vertex_to_cell_w_inv_distance();

      const cs_weight_t *w = _weights_vtc[CS_VERTEX_TO_CELL_SHEPARD];

      cs_real_t *c_w = nullptr;
      if (v_weight != nullptr) {
        BFT_MALLOC(c_w, n_cells, cs_real_t);
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          c_w[c_id] = 0;
      }

      if (v_weight == nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              c_var[c_id*var_dim + k] += v_var[v_id*var_dim + k] * w[j];
          }
        }
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              c_var[c_id*var_dim + k] +=   v_var[v_id*var_dim + k] * w[j]
                                         * v_weight[v_id];
            c_w[c_id] += w[j] * v_weight[v_id];
          }
        }
      }

      if (v_weight != nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          for (cs_lnum_t k = 0; k < var_dim; k++)
            c_var[c_id*var_dim+k] /= c_w[c_id];
        BFT_FREE(c_w);
      }

    }
    break;

  case CS_VERTEX_TO_CELL_LR:
    {
      if (! _set_vtc[CS_VERTEX_TO_CELL_LR])
        _vertex_to_cell_f_lsq();

      cs_real_t  *rhs;
      cs_lnum_t  rhs_size = n_cells*4*var_dim;
      BFT_MALLOC(rhs, rhs_size, cs_real_t);
      for (cs_lnum_t i = 0; i < rhs_size; i++)
        rhs[i] = 0;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t *c_coo = mq->cell_cen + c_id*3;
        cs_lnum_t s_id = c2v_idx[c_id];
        cs_lnum_t e_id = c2v_idx[c_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t v_id = c2v_ids[j];
          const cs_real_t *_v_var = v_var + v_id*var_dim;
          const cs_real_t *v_coo = m->vtx_coord + v_id*3;
          cs_real_t r_coo[3]
            = {v_coo[0]-c_coo[0], v_coo[1]-c_coo[1], v_coo[2]-c_coo[2]};
          for (cs_lnum_t k = 0; k < var_dim; k++) {
            cs_real_t  *_rhs = rhs + c_id*4*var_dim + k*4;
            _rhs[0] += r_coo[0] * _v_var[k];
            _rhs[1] += r_coo[1] * _v_var[k];
            _rhs[2] += r_coo[2] * _v_var[k];
            _rhs[3] += _v_var[k];
          }
        }
      }

      const cs_weight_t *ldlt = _weights_vtc[CS_VERTEX_TO_CELL_LR];

#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t  *_ldlt = ldlt + c_id*10;
        for (cs_lnum_t k = 0; k < var_dim; k++) {
          const cs_real_t  *_rhs = rhs + c_id*4*var_dim + k*4;
          c_var[c_id*var_dim + k] = cs_math_sym_44_partial_solve_ldlt(_ldlt, _rhs);
        }
      }

      BFT_FREE(rhs);
    }
    break;
  default:
    break;
  }
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
cs_vertex_to_cell_free(void)
{
  for (int i = 0; i < 3; i++)
    BFT_FREE(_weights_vtc[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate vertex values to cell values.
 *
 * \param[in]       method            interpolation method
 * \param[in]       verbosity         verbosity level
 * \param[in]       var_dim           variable dimension
 * \param[in]       ignore_rot_perio  if true, ignore periodicity of rotation
 * \param[in]       v_weight          vertex weight, or nullptr
 * \param[in]       v_var             base vertex-based variable
 * \param[out]      c_var             cell-based variable
 */
/*----------------------------------------------------------------------------*/

void
cs_vertex_to_cell(cs_vertex_to_cell_type_t   method,
                  int                        verbosity,
                  cs_lnum_t                  var_dim,
                  const cs_real_t *restrict  v_weight,
                  const cs_real_t *restrict  v_var,
                  cs_real_t *restrict        c_var)
{
  CS_UNUSED(verbosity);

  if (var_dim == 1)
    _vertex_to_cell_scalar(method,
                           verbosity,
                           v_weight,
                           v_var,
                           c_var);

  else
    _vertex_to_cell_strided(method,
                            verbosity,
                            var_dim,
                            v_weight,
                            v_var,
                            c_var);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
