/*============================================================================
 * Gradient reconstruction, CUDA implementations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Initialize RHS with null values
 *----------------------------------------------------------------------------*/

__global__ static void
_init_rhs_v2(cs_lnum_t         n_cells_g,
           cs_real_33_t      *restrict _rhs)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells_g)
    return;

  rhs[c_id] = 0.0;
}

__global__ static void
_compute_rhs_lsq_v_i_face_v2(cs_lnum_t            n_i_faces,
                          const cs_lnum_t      *restrict i_face_cells,
                          const cs_real_t    *restrict cell_f_cen,
                          cs_real_33_t         *restrict _rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t         *restrict weight,
                          const cs_real_t      *restrict c_weight)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _weight1, _weight2, _denom, _pond, pfac;
  cs_lnum_t c_id1, c_id2;

  c_id1 = i_face_cells[f_id*2];
  c_id2 = i_face_cells[f_id*2 + 1];

  dc[0] = cell_f_cen[c_id2*3] - cell_f_cen[c_id1*3];
  dc[1] = cell_f_cen[c_id2*3 + 1] - cell_f_cen[c_id1*3 + 1];
  dc[2] = cell_f_cen[c_id2*3 + 2] - cell_f_cen[c_id1*3 + 2];

  ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
  
  if (c_weight == NULL){
    _weight1 = 1.;
    _weight2 = 1.;
  }
  else{
    _pond = weight[f_id];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
    _weight1 = c_weight[c_id1] * _denom;
    _weight2 = c_weight[c_id2] * _denom;
  }
  
  for(cs_lnum_t i = 0; i < 3; i++){
    pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
    for(cs_lnum_t j = 0; j < 3; j++){
      fctb[j] = dc[j] * pfac;
      atomicAdd(&rhs[c_id1*3*3 + i*3 + j], _weight2 * fctb[j]);
      atomicAdd(&rhs[c_id2*3*3 + i*3 + j], _weight1 * fctb[j]);
    }
  }
}

__global__ static void
_compute_rhs_lsq_v_i_face_v2cf(cs_lnum_t            size,
                          const cs_lnum_2_t      *restrict _i_face_cells,
                          const cs_real_3_t    *restrict _cell_f_cen,
                          cs_real_33_t         *restrict _rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t         *restrict weight,
                          const cs_real_t      *restrict c_weight)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_lnum_t *i_face_cells = (cs_lnum_t *) _i_face_cells;
  cs_real_t *cell_f_cen = (cs_real_t *) _cell_f_cen;
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= size){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _weight1, _weight2, _denom, _pond, pfac;
  cs_lnum_t c_id1, c_id2;

  c_id1 = i_face_cells[f_id*2];
  c_id2 = i_face_cells[f_id*2 + 1];

  dc[0] = cell_f_cen[c_id2*3] - cell_f_cen[c_id1*3];
  dc[1] = cell_f_cen[c_id2*3 + 1] - cell_f_cen[c_id1*3 + 1];
  dc[2] = cell_f_cen[c_id2*3 + 2] - cell_f_cen[c_id1*3 + 2];

  ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
  
  if (c_weight == NULL){
    _weight1 = 1.;
    _weight2 = 1.;
  }
  else{
    _pond = weight[f_id];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
    _weight1 = c_weight[c_id1] * _denom;
    _weight2 = c_weight[c_id2] * _denom;
  }

  using Cell = AtomicCell<cs_real_t, 3, 3>;
  Cell _rhs1, _rhs2;
  
  for(cs_lnum_t i = 0; i < 3; i++){
    pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
    for(cs_lnum_t j = 0; j < 3; j++){
      fctb[j] = dc[j] * pfac;
      //_rhs1[i][j].get() += _weight2 * fctb[j];
      //_rhs2[i][j].get() += _weight1 * fctb[j];
      atomicAdd(&rhs[c_id1*3*3 + i*3 + j], _weight2 * fctb[j]);
      atomicAdd(&rhs[c_id2*3*3 + i*3 + j], _weight1 * fctb[j]);
    }
  }
  //reinterpret_cast<Cell*>(&rhs[c_id1*3*3][0][0])->atomic_add(_rhs1);
  //reinterpret_cast<Cell*>(&rhs[c_id2*3*3][0][0])->atomic_add(_rhs2);
}

__global__ static void
_compute_rhs_lsq_v_b_face_v2(cs_lnum_t           n_b_faces,
                            const cs_lnum_t      *restrict b_face_cells,
                            const cs_real_3_t    *restrict cell_f_cen,
                            const cs_real_3_t    *restrict b_face_normal,
                            cs_real_33_t            *restrict _rhs,
                            const cs_real_3_t    *restrict pvar,
                            const cs_real_t      *restrict b_dist,
                            const cs_real_33_t   *restrict coefbv,
                            const cs_real_3_t    *restrict coefav,
                            const int            inc)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_b_faces){
    return;
  }

  cs_lnum_t c_id1;
  cs_real_t n_d_dist[3], d_b_dist, pfac, norm, inverse_norm;

  c_id1 = b_face_cells[f_id];

  cs_math_3_normalize_cuda(b_face_normal[f_id], n_d_dist);

  d_b_dist = 1. / b_dist[f_id];

  /* Normal divided by b_dist */
  n_d_dist[0] *= d_b_dist;
  n_d_dist[1] *= d_b_dist;
  n_d_dist[2] *= d_b_dist;

  for (cs_lnum_t i = 0; i < 3; i++) {
    pfac =   coefav[f_id][i]*inc
          + ( coefbv[f_id][0][i] * pvar[c_id1][0]
            + coefbv[f_id][1][i] * pvar[c_id1][1]
            + coefbv[f_id][2][i] * pvar[c_id1][2]
            - pvar[c_id1][i]);

    atomicAdd(&rhs[c_id1*3*3 + i*3], n_d_dist[0] * pfac);
    atomicAdd(&rhs[c_id1*3*3 + i*3 + 1], n_d_dist[1] * pfac);
    atomicAdd(&rhs[c_id1*3*3 + i*3 + 2], n_d_dist[2] * pfac); 
  }
}

__global__ static void
_compute_gradient_lsq_v_v2(cs_lnum_t           n_cells_g,
                        cs_real_33_t        *restrict _gradv,
                        cs_real_33_t        *restrict _rhs,
                        cs_cocg_6_t         *restrict cocg)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_real_t *gradv = (cs_real_t *) _gradv;
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells_g) 
    return;

  for(cs_lnum_t i = 0; i < 3; i++){
    gradv[c_id*3*3 + i*3] =   rhs[c_id*3*3 + i*3] * cocg[c_id][0]
                          + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][3]
                          + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][5];

    gradv[c_id*3*3 + i*3 + 1] =   rhs[c_id*3*3 + i*3] * cocg[c_id][3]
                        + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][1]
                        + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][4];

    gradv[c_id*3*3 + i*3 + 2] =   rhs[c_id*3*3 + i*3] * cocg[c_id][5]
                        + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][4]
                        + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][2];
  }
}

__global__ static void
_compute_gradient_lsq_v_v4(cs_lnum_t           n_cells,
                        cs_real_33_t        *restrict gradv_m,
                        cs_real_33_t        *restrict rhs_m,
                        cs_cocg_6_t         *restrict cocg)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells) 
    return;

  cs_real_t *rhs = (cs_real_t *) rhs_m;
  cs_real_t *gradv = (cs_real_t *) gradv_m;

  for(cs_lnum_t i = 0; i < 3; i++){
    gradv[c_id*3*3 + i*3] =   rhs[c_id*3*3 + i*3] * cocg[c_id][0]
                          + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][3]
                          + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][5];

    gradv[c_id*3*3 + i*3 + 1] =   rhs[c_id*3*3 + i*3] * cocg[c_id][3]
                        + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][1]
                        + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][4];

    gradv[c_id*3*3 + i*3 + 2] =   rhs[c_id*3*3 + i*3] * cocg[c_id][5]
                        + rhs[c_id*3*3 + i*3 + 1] * cocg[c_id][4]
                        + rhs[c_id*3*3 + i*3 + 2] * cocg[c_id][2];
  }
}
