/*============================================================================
 * Cell to vertex interpolation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_cell_to_vertex.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cell_to_vertex.c
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

bool          _set[3] = {false, false, false};
cs_weight_t  *_weights[3][2] = {{NULL, NULL}, {NULL, NULL}, {NULL, NULL}};

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
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (m00, m10, m11, m20, m21, m22, m30, m31, m32, m33) in
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33) out
 */
/*----------------------------------------------------------------------------*/

static inline void
_sym_44_factor_ldlt(cs_real_t  ldlt[10])
{
  /* Factorization */

  // j=0
  const cs_real_t  d00 = ldlt[0]; // m00
  assert(fabs(d00) > cs_math_zero_threshold);

  const cs_real_t  f00 = 1. / d00;
  const cs_real_t  l10 = ldlt[1] * f00; // m01
  const cs_real_t  l20 = ldlt[3] * f00; // m02
  const cs_real_t  l30 = ldlt[6] * f00; // m03

  // j=1
  const cs_real_t  d11 = ldlt[2] - l10*l10*d00; // m11
  assert(fabs(d11) > cs_math_zero_threshold);
  const cs_real_t  f11 = 1. / d11;
  const cs_real_t  l21 = (ldlt[4] - l20*d00*l10) * f11; // m12
  const cs_real_t  l31 = (ldlt[7] - l30*d00*l10) * f11; // m13

  // j=2
  const cs_real_t  d22 = ldlt[5] - l20*d00*l20 - l21*d11*l21;  // m22
  assert(fabs(d22) > cs_math_zero_threshold);
  const cs_real_t  f22 = 1. / d22;
  const cs_real_t  l32 = (ldlt[8] - l30*d00*l20 - l31*d11*l21) * f22;  // m23

  // j=3
  const cs_real_t  d33 = ldlt[9] - l30*d00*l30 - l31*d11*l31 - l32*d22*l32; // m33
  assert(fabs(d33) > cs_math_zero_threshold);
  const cs_real_t  f33 = 1. / d33;

  ldlt[0] = f00;
  ldlt[1] = l10;
  ldlt[2] = f11;
  ldlt[3] = l20;
  ldlt[4] = l21;
  ldlt[5] = f22;
  ldlt[6] = l30;
  ldlt[7] = l31;
  ldlt[8] = l32;
  ldlt[9] = f33;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33)
 * \param[in]       rhs   right-hand side
 * \param[out]      x     solution
 */
/*----------------------------------------------------------------------------*/

static void
_sym_44_solve_ldlt(const cs_real_t  ldlt[10],
                   const cs_real_t  rhs[4],
                   cs_real_t        x[4])
{
  /* f00, f11, f22, f33, l32, l31, l30, l21, l20, l10
     0    1    2    3    4    5    6    7    8    9   */

  x[0] = rhs[0];
  x[1] = rhs[1] - x[0]*ldlt[1];
  x[2] = rhs[2] - x[0]*ldlt[3] - x[1]*ldlt[4];
  x[3] = rhs[3] - x[0]*ldlt[6] - x[1]*ldlt[7] - x[2]*ldlt[8];

  x[3] = x[3]*ldlt[9];
  x[2] = x[2]*ldlt[5] - ldlt[8]*x[3];
  x[1] = x[1]*ldlt[2] - ldlt[7]*x[3] - ldlt[4]*x[2];
  x[0] = x[0]*ldlt[0] - ldlt[6]*x[3] - ldlt[3]*x[2] - ldlt[1]*x[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute multipliers (pseudo-weights) for unweighted case
 *
 * In this case the weights directly match the number of adjacent cells.
 *
 * \param[in]  ignore_r_tr  ignore periodicity with rotation if true
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_w_unweighted(bool  ignore_r_tr)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0];
  int *w_sum;
  BFT_REALLOC(w, n_vertices, cs_weight_t);
  BFT_MALLOC(w_sum, n_vertices, int);

  _set[CS_CELL_TO_VERTEX_UNWEIGHTED] = true;
  _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0] = w;

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    w_sum[v_id] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      w_sum[v_id] += 1;
    }

  }

  if (m->vtx_interfaces != NULL)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            1,
                            true,
                            CS_INT_TYPE,
                            ignore_r_tr,
                            w_sum);

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    w[v_id] = 1. / w_sum[v_id];

  BFT_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute weights based on inversed distance
 *
 * \param[in]  ignore_r_tr  ignore periodicity with rotation if true
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_w_inv_distance(bool  ignore_r_tr)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_SHEPARD][0];
  cs_weight_t *wb = _weights[CS_CELL_TO_VERTEX_SHEPARD][1];
  cs_real_t   *w_sum;
  BFT_REALLOC(w, c2v_idx[n_cells], cs_weight_t);
  BFT_REALLOC(wb, f2v_idx[n_b_faces], cs_weight_t);
  BFT_MALLOC(w_sum, n_vertices, cs_real_t);

  _set[CS_CELL_TO_VERTEX_SHEPARD] = true;
  _weights[CS_CELL_TO_VERTEX_SHEPARD][0] = w;
  _weights[CS_CELL_TO_VERTEX_SHEPARD][1] = wb;

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    w_sum[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t *c_coo = mq->cell_cen + c_id*3;

    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      cs_real_t *v_coo = m->vtx_coord + v_id*3;
      cs_real_t d = cs_math_3_distance(c_coo, v_coo);
      if (d <= DBL_MIN) {
        w[j] = 1;
        w_sum[v_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1./d;
        w[j] = _w;
        w_sum[v_id] += _w;
      }
    }

  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_real_t *f_coo = mq->b_face_cog + f_id*3;

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = f2v_ids[j];
      cs_real_t *v_coo = m->vtx_coord + v_id*3;
      cs_real_t d = cs_math_3_distance(f_coo, v_coo);
      if (d <= DBL_MIN) {
        wb[j] = 1;
        w_sum[v_id] = HUGE_VAL;
      }
      else {
        cs_real_t _w = 1./d;
        wb[j] = _w;
        w_sum[v_id] += _w;
      }
    }

  }

  if (m->vtx_interfaces != NULL)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            1,
                            true,
                            CS_REAL_TYPE,
                            ignore_r_tr,
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

  BFT_FREE(w_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute factorization based on least squares
 *
 * \param[in]  ignore_r_tr  ignore periodicity with rotation if true
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_f_lsq(bool  ignore_r_tr)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_glob_mesh_adjacencies->c2v;

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  cs_lnum_t  w_size = n_vertices*10;

  cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_LR][0];
  BFT_MALLOC(w, w_size, cs_real_t);

  _set[CS_CELL_TO_VERTEX_LR] = true;
  _weights[CS_CELL_TO_VERTEX_LR][0] = w;

# pragma omp parallel for if(w_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < w_size; i++)
    w[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t *c_coo = mq->cell_cen + c_id*3;
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      cs_real_t  *_a = w + v_id*10;
      _a[0] += c_coo[0] * c_coo[0]; // a00
      _a[1] += c_coo[1] * c_coo[0]; // a10
      _a[2] += c_coo[1] * c_coo[1]; // a11
      _a[3] += c_coo[2] * c_coo[0]; // a20
      _a[4] += c_coo[2] * c_coo[1]; // a21
      _a[5] += c_coo[2] * c_coo[2]; // a22
      _a[6] += c_coo[0];            // a30
      _a[7] += c_coo[1];            // a31
      _a[8] += c_coo[2];            // a32
      _a[9] += 1;                   // a33
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    const cs_real_t *f_coo = mq->b_face_cog + f_id*3;
    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = f2v_ids[j];
      cs_real_t  *_a = w + v_id*10;
      _a[0] += f_coo[0] * f_coo[0]; // a00
      _a[1] += f_coo[1] * f_coo[0]; // a10
      _a[2] += f_coo[1] * f_coo[1]; // a11
      _a[3] += f_coo[2] * f_coo[0]; // a20
      _a[4] += f_coo[2] * f_coo[1]; // a21
      _a[5] += f_coo[2] * f_coo[2]; // a22
      _a[6] += f_coo[0];            // a30
      _a[7] += f_coo[1];            // a31
      _a[8] += f_coo[2];            // a32
      _a[9] += 1;                   // a33
    }
  }

  if (m->vtx_interfaces != NULL)
    cs_interface_set_sum_tr(m->vtx_interfaces,
                            m->n_vertices,
                            10,
                            true,
                            CS_REAL_TYPE,
                            ignore_r_tr,
                            w);

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    _sym_44_factor_ldlt(w + v_id*10);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate cell values to vertex values for a scalar arrray.
 *
 * \param[in]  method       interpolation method
 * \param[in]  verbosity    verbosity level
 * \param[in]  ignore_r_tr  ignore periodicity with rotation if true
 * \param[in]  c_weight     cell weight, or NULL
 * \param[in]  c_var        base cell-based variable
 * \param[in]  b_var        base boundary-face values, or NULL
 * \param[out] v_var        vertex-based variable
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_scalar(cs_cell_to_vertex_type_t   method,
                       int                        verbosity,
                       bool                       ignore_r_tr,
                       const cs_real_t            c_weight[restrict],
                       const cs_real_t            c_var[restrict],
                       const cs_real_t            b_var[restrict],
                       cs_real_t                  v_var[restrict])
{
  CS_UNUSED(verbosity);

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_mesh_adjacencies_cell_vertices();

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

# pragma omp parallel for if(n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    v_var[v_id] = 0;

  switch(method) {

  case CS_CELL_TO_VERTEX_UNWEIGHTED:
    {
      if (c_weight == NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            v_var[v_id] += c_var[c_id];
          }
        }

        if (m->vtx_interfaces != NULL)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  m->n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_var);

        if (! _set[CS_CELL_TO_VERTEX_UNWEIGHTED])
          _cell_to_vertex_w_unweighted(ignore_r_tr);

        const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0];
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_var[v_id] *= w[v_id];
      }
      else {
        cs_real_t *v_w;
        BFT_MALLOC(v_w, n_vertices, cs_real_t);
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_w[v_id] = 0;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            v_var[v_id] += c_var[c_id]*c_weight[c_id];
            v_w[v_id] += c_weight[c_id];
          }
        }

        if (m->vtx_interfaces != NULL) {
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_var);
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_w);
        }

        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_var[v_id] /= v_w[v_id];

        BFT_FREE(v_w);
      }
    }
    break;

  case CS_CELL_TO_VERTEX_SHEPARD:
    {
      if (! _set[CS_CELL_TO_VERTEX_SHEPARD])
        _cell_to_vertex_w_inv_distance(ignore_r_tr);

      const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_SHEPARD][0];

      cs_real_t *v_w = NULL;
      if (c_weight != NULL) {
        BFT_MALLOC(v_w, n_vertices, cs_real_t);
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_w[v_id] = 0;
      }

      if (c_weight == NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            v_var[v_id] += c_var[c_id] * w[j];
          }
        }
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            v_var[v_id] += c_var[c_id] * w[j] * c_weight[c_id];
            v_w[v_id] += c_weight[c_id];
          }
        }
      }

      const cs_weight_t *wb = _weights[CS_CELL_TO_VERTEX_SHEPARD][1];

      if (c_weight == NULL) {

        if (b_var == NULL) {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            cs_real_t _b_var = c_var[c_id];
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              v_var[v_id] += _b_var * wb[j];
            }
          }
        }
        else {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
            v_var[v_id] += b_var[f_id] * wb[j];
            }
          }
        }

      }
      else { /* c_weight != NULL */

        if (b_var == NULL) {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            cs_real_t _b_var = c_var[c_id];
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              v_var[v_id] += _b_var * wb[j] * c_weight[c_id];
              v_w[v_id] += c_weight[c_id];
            }
          }
        }
        else {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              v_var[v_id] += b_var[f_id] * wb[j] * c_weight[c_id];
              v_w[v_id] += c_weight[c_id];
            }
          }
        }
      }

      if (m->vtx_interfaces != NULL)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                1,
                                true,
                                CS_REAL_TYPE,
                                ignore_r_tr,
                                v_var);

      if (c_weight != NULL) {
        if (m->vtx_interfaces != NULL)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_w);
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_var[v_id] /= v_w[v_id];
        BFT_FREE(v_w);
      }

    }
    break;

  case CS_CELL_TO_VERTEX_LR:
    {
      if (! _set[CS_CELL_TO_VERTEX_LR])
        _cell_to_vertex_f_lsq(ignore_r_tr);

      cs_real_t  *rhs;
      cs_lnum_t  rhs_size = n_vertices*4;
      BFT_MALLOC(rhs, rhs_size, cs_real_t);
      for (cs_lnum_t i = 0; i < rhs_size; i++)
        rhs[i] = 0;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t *c_coo = mq->cell_cen + c_id*3;
        cs_real_t _c_var = c_var[c_id];
        cs_lnum_t s_id = c2v_idx[c_id];
        cs_lnum_t e_id = c2v_idx[c_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t v_id = c2v_ids[j];
          cs_real_t  *_rhs = rhs + v_id*4;
          _rhs[0] += c_coo[0] * _c_var;
          _rhs[1] += c_coo[1] * _c_var;
          _rhs[2] += c_coo[2] * _c_var;
          _rhs[3] += c_var[c_id];
        }
      }

      if (b_var == NULL) {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          const cs_real_t *f_coo = mq->b_face_cog + f_id*3;
          cs_lnum_t c_id = m->b_face_cells[f_id];
          cs_real_t _b_var = c_var[c_id];
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            cs_real_t  *_rhs = rhs + v_id*4;
            _rhs[0] += f_coo[0] * _b_var;
            _rhs[1] += f_coo[1] * _b_var;
            _rhs[2] += f_coo[2] * _b_var;
            _rhs[3] += _b_var;
          }
        }
      }
      else {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          const cs_real_t *f_coo = mq->b_face_cog + f_id*3;
          cs_real_t _b_var = b_var[f_id];
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            cs_real_t  *_rhs = rhs + v_id*4;
            _rhs[0] += f_coo[0] * _b_var;
            _rhs[1] += f_coo[1] * _b_var;
            _rhs[2] += f_coo[2] * _b_var;
            _rhs[3] += _b_var;
          }
        }
      }

      if (m->vtx_interfaces != NULL)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                4,
                                true,
                                CS_REAL_TYPE,
                                ignore_r_tr,
                                rhs);

      const cs_weight_t *ldlt = _weights[CS_CELL_TO_VERTEX_LR][0];

#     pragma omp parallel for if(n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
        const cs_real_t *v_coo = m->vtx_coord + v_id*3;
        const cs_real_t  *_ldlt = ldlt + v_id*10;
        const cs_real_t  *_rhs = rhs + v_id*4;
        cs_real_t x[4];
        _sym_44_solve_ldlt(_ldlt, _rhs, x);
        v_var[v_id] = (  x[0]*v_coo[0]
                       + x[1]*v_coo[1]
                       + x[2]*v_coo[2]
                       + x[3]);
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
 * \brief  Interpolate cell values to vertex values for a strided arrray.
 *
 * \param[in]   method      interpolation method
 * \param[in]   verbosity   verbosity level
 * \param[in]   var_dim     varible dimension
 * \param[in]   ignore_r_tr  ignore periodicity with rotation if true
 * \param[in]   c_weight    cell weight, or NULL
 * \param[in]   c_var       base cell-based variable
 * \param[in]   b_var       base boundary-face values, or NULL
 * \param[out]  v_var       vertex-based variable
 */
/*----------------------------------------------------------------------------*/

static void
_cell_to_vertex_strided(cs_cell_to_vertex_type_t   method,
                        int                        verbosity,
                        cs_lnum_t                  var_dim,
                        bool                       ignore_r_tr,
                        const cs_real_t            c_weight[restrict],
                        const cs_real_t            c_var[restrict],
                        const cs_real_t            b_var[restrict],
                        cs_real_t                  v_var[restrict])
{
  CS_UNUSED(verbosity);

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_adjacency_t  *c2v = cs_mesh_adjacencies_cell_vertices();

  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

  const cs_lnum_t *f2v_idx = m->b_face_vtx_idx;
  const cs_lnum_t *f2v_ids = m->b_face_vtx_lst;

  const cs_lnum_t n_v_values = n_vertices*var_dim;

# pragma omp parallel for if(n_v_values > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < n_v_values; v_id++)
    v_var[v_id] = 0;

  switch(method) {

  case CS_CELL_TO_VERTEX_UNWEIGHTED:
    {
      if (c_weight == NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              v_var[v_id*var_dim + k] += c_var[c_id*var_dim + k];
          }
        }

        if (m->vtx_interfaces != NULL)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  m->n_vertices,
                                  var_dim,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_var);

        if (! _set[CS_CELL_TO_VERTEX_UNWEIGHTED])
          _cell_to_vertex_w_unweighted(ignore_r_tr);

        const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_UNWEIGHTED][0];
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
          for (cs_lnum_t k = 0; k < var_dim; k++)
            v_var[v_id*var_dim + k] *= w[v_id];
        }
      }
      else {
        cs_real_t *v_w;
        BFT_MALLOC(v_w, n_vertices, cs_real_t);
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
          v_w[v_id] = 0;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              v_var[v_id*var_dim + k] += c_var[c_id*var_dim + k] * c_weight[c_id];
            v_w[v_id] += c_weight[c_id];
          }
        }

        if (m->vtx_interfaces != NULL) {
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  var_dim,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_var);
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_w);
        }
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
          for (cs_lnum_t k = 0; k < var_dim; k++)
            v_var[v_id*var_dim + k] /= v_w[v_id];
        }

        BFT_FREE(v_w);
      }
    }
    break;

  case CS_CELL_TO_VERTEX_SHEPARD:
    {
      if (! _set[CS_CELL_TO_VERTEX_SHEPARD])
        _cell_to_vertex_w_inv_distance(ignore_r_tr);

      const cs_weight_t *w = _weights[CS_CELL_TO_VERTEX_SHEPARD][0];

      cs_real_t *v_w = NULL;
      if (c_weight != NULL) {
        BFT_MALLOC(v_w, n_vertices, cs_real_t);
        for (cs_lnum_t v_id = 0; v_id < n_v_values; v_id++)
          v_w[v_id] = 0;
      }

      if (c_weight == NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_lnum_t s_id = c2v_idx[c_id];
          cs_lnum_t e_id = c2v_idx[c_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = c2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++)
              v_var[v_id*var_dim + k] += c_var[c_id*var_dim + k] * w[j];
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
              v_var[v_id*var_dim + k] +=   c_var[c_id*var_dim + k] * w[j]
                                         * c_weight[c_id];
            v_w[v_id] += c_weight[c_id];
          }
        }
      }

      const cs_weight_t *wb = _weights[CS_CELL_TO_VERTEX_SHEPARD][1];

      if (c_weight == NULL) {

        if (b_var == NULL) {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            const cs_real_t *_b_var = c_var + c_id*var_dim;
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < var_dim; k++)
                v_var[v_id*var_dim + k] += _b_var[k] * wb[j];
            }
          }
        }
        else {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < var_dim; k++)
                v_var[v_id*var_dim+k] += b_var[f_id*var_dim+k] * wb[j];
            }
          }
        }

      }
      else { /* c_weight != NULL */

        if (b_var == NULL) {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            const cs_real_t *_b_var = c_var + c_id*var_dim;
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < var_dim; k++)
                v_var[v_id*var_dim + k] += _b_var[k] * wb[j] * c_weight[c_id];
              v_w[v_id] += c_weight[c_id];
            }
          }
        }
        else {
          for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
            cs_lnum_t c_id = m->b_face_cells[f_id];
            cs_lnum_t s_id = f2v_idx[f_id];
            cs_lnum_t e_id = f2v_idx[f_id+1];
            for (cs_lnum_t j = s_id; j < e_id; j++) {
              cs_lnum_t v_id = f2v_ids[j];
              for (cs_lnum_t k = 0; k < var_dim; k++)
                v_var[v_id*var_dim + k] +=   b_var[f_id*var_dim + k] * wb[j]
                                           * c_weight[c_id];
              v_w[v_id] += c_weight[c_id];
            }
          }
        }
      }

      if (m->vtx_interfaces != NULL)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                var_dim,
                                true,
                                CS_REAL_TYPE,
                                ignore_r_tr,
                                v_var);

      if (c_weight != NULL) {
        if (m->vtx_interfaces != NULL)
          cs_interface_set_sum_tr(m->vtx_interfaces,
                                  n_vertices,
                                  1,
                                  true,
                                  CS_REAL_TYPE,
                                  ignore_r_tr,
                                  v_w);
        for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
          for (cs_lnum_t k = 0; k < var_dim; k++)
            v_var[v_id*var_dim+k] /= v_w[v_id];
        }
        BFT_FREE(v_w);
      }

    }
    break;

  case CS_CELL_TO_VERTEX_LR:
    {
      if (! _set[CS_CELL_TO_VERTEX_LR])
        _cell_to_vertex_f_lsq(ignore_r_tr);

      cs_real_t  *rhs;
      cs_lnum_t  rhs_size = n_vertices*4*var_dim;
      BFT_MALLOC(rhs, rhs_size, cs_real_t);
      for (cs_lnum_t i = 0; i < rhs_size; i++)
        rhs[i] = 0;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t *c_coo = mq->cell_cen + c_id*3;
        const cs_real_t *_c_var = c_var + c_id*var_dim;
        cs_lnum_t s_id = c2v_idx[c_id];
        cs_lnum_t e_id = c2v_idx[c_id+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_lnum_t v_id = c2v_ids[j];
          for (cs_lnum_t k = 0; k < var_dim; k++) {
            cs_real_t  *_rhs = rhs + v_id*4*var_dim + k*4;
            _rhs[0] += c_coo[0] * _c_var[k];
            _rhs[1] += c_coo[1] * _c_var[k];
            _rhs[2] += c_coo[2] * _c_var[k];
            _rhs[3] += _c_var[k];
          }
        }
      }

      if (b_var == NULL) {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          const cs_real_t *f_coo = mq->b_face_cog + f_id*3;
          cs_lnum_t c_id = m->b_face_cells[f_id];
          const cs_real_t *_b_var = c_var + c_id*var_dim;
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++) {
              cs_real_t  *_rhs = rhs + v_id*4*var_dim + k*4;
              _rhs[0] += f_coo[0] * _b_var[k];
              _rhs[1] += f_coo[1] * _b_var[k];
              _rhs[2] += f_coo[2] * _b_var[k];
              _rhs[3] += _b_var[k];
            }
          }
        }
      }
      else {
        for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
          const cs_real_t *f_coo = mq->b_face_cog + f_id*3;
          const cs_real_t *_b_var = b_var + f_id*var_dim;
          cs_lnum_t s_id = f2v_idx[f_id];
          cs_lnum_t e_id = f2v_idx[f_id+1];
          for (cs_lnum_t j = s_id; j < e_id; j++) {
            cs_lnum_t v_id = f2v_ids[j];
            for (cs_lnum_t k = 0; k < var_dim; k++) {
              cs_real_t  *_rhs = rhs + v_id*4*var_dim + k*4;
              _rhs[0] += f_coo[0] * _b_var[k];
              _rhs[1] += f_coo[1] * _b_var[k];
              _rhs[2] += f_coo[2] * _b_var[k];
              _rhs[3] += _b_var[k];
            }
          }
        }
      }

      if (m->vtx_interfaces != NULL)
        cs_interface_set_sum_tr(m->vtx_interfaces,
                                m->n_vertices,
                                4*var_dim,
                                true,
                                CS_REAL_TYPE,
                                ignore_r_tr,
                                rhs);

      const cs_weight_t *ldlt = _weights[CS_CELL_TO_VERTEX_LR][0];

#     pragma omp parallel for if(n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
        const cs_real_t *v_coo = m->vtx_coord + v_id*3;
        const cs_real_t  *_ldlt = ldlt + v_id*10;
        for (cs_lnum_t k = 0; k < var_dim; k++) {
          const cs_real_t  *_rhs = rhs + v_id*4*var_dim + k*4;
          cs_real_t x[4];
          _sym_44_solve_ldlt(_ldlt, _rhs, x);
          v_var[v_id*var_dim + k] = (  x[0]*v_coo[0]
                                     + x[1]*v_coo[1]
                                     + x[2]*v_coo[2]
                                     + x[3]);
        }
      }

      BFT_FREE(rhs);
    }
    break;
  default:
    break;
  }
}

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
    BFT_FREE(_weights[i][j]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate cell values to vertex values.
 *
 * \param[in]       method            interpolation method
 * \param[in]       verbosity         verbosity level
 * \param[in]       var_dim           varible dimension
 * \param[in]       ignore_rot_perio  if true, ignore periodicity of rotation
 * \param[in]       c_weight          cell weight, or NULL
 * \param[in]       c_var             base cell-based variable
 * \param[in]       b_var             base boundary-face values, or NULL
 * \param[out]      v_var             vertex-based variable
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_to_vertex(cs_cell_to_vertex_type_t   method,
                  int                        verbosity,
                  cs_lnum_t                  var_dim,
                  bool                       ignore_rot_perio,
                  const cs_real_t            c_weight[restrict],
                  const cs_real_t            c_var[restrict],
                  const cs_real_t            b_var[restrict],
                  cs_real_t                  v_var[restrict])
{
  CS_UNUSED(verbosity);

  if (var_dim == 1)
    _cell_to_vertex_scalar(method,
                           verbosity,
                           ignore_rot_perio,
                           c_weight,
                           c_var,
                           b_var,
                           v_var);

  else
    _cell_to_vertex_strided(method,
                            verbosity,
                            var_dim,
                            ignore_rot_perio,
                            c_weight,
                            c_var,
                            b_var,
                            v_var);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
