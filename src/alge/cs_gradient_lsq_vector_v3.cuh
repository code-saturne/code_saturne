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


__global__ static void
_compute_rhs_lsq_v_i_face_v3(cs_lnum_t            n_i_faces,
                          const cs_lnum_2_t      *restrict i_face_cells,
                          const cs_real_3_t    *restrict cell_f_cen,
                          cs_real_33_t         *restrict rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t         *restrict weight,
                          const cs_real_t      *restrict c_weight)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _weight1, _weight2, _denom, _pond, pfac;
  cs_lnum_t c_id1, c_id2;

  size_t f_id1 = f_id / (3*3);
  size_t i = (f_id / 3) % 3;
  size_t j = f_id % 3;

  c_id1 = i_face_cells[f_id1][0];
  c_id2 = i_face_cells[f_id1][1];

  dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
  dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
  dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

  ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
  
  if (c_weight == NULL){
    _weight1 = 1.;
    _weight2 = 1.;
  }
  else{
    _pond = weight[f_id1];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
    _weight1 = c_weight[c_id1] * _denom;
    _weight2 = c_weight[c_id2] * _denom;
  }
  
  //for(cs_lnum_t i = 0; i < 3; i++){
    pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
    //for(cs_lnum_t j = 0; j < 3; j++){
      fctb[j] = dc[j] * pfac;
      atomicAdd(&rhs[c_id1][i][j], _weight2 * fctb[j]);
      atomicAdd(&rhs[c_id2][i][j], _weight1 * fctb[j]);
    //}
  //}
}

__global__ static void
_compute_rhs_lsq_v_i_face_v3cf(cs_lnum_t            size,
                          const cs_lnum_2_t      *restrict i_face_cells,
                          const cs_real_3_t    *restrict cell_f_cen,
                          cs_real_33_t         *restrict rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t         *restrict weight,
                          const cs_real_t      *restrict c_weight)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= size){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _weight1, _weight2, _denom, _pond, pfac;
  cs_lnum_t c_id1, c_id2;

  size_t f_id1 = f_id / (3*3);
  size_t i = (f_id / 3) % 3;
  size_t j = f_id % 3;

  c_id1 = i_face_cells[f_id1][0];
  c_id2 = i_face_cells[f_id1][1];

  dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
  dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
  dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

  ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
  
  if (c_weight == NULL){
    _weight1 = 1.;
    _weight2 = 1.;
  }
  else{
    _pond = weight[f_id1];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
    _weight1 = c_weight[c_id1] * _denom;
    _weight2 = c_weight[c_id2] * _denom;
  }

  using Cell = AtomicCell<cs_real_t>;
  
  //for(cs_lnum_t i = 0; i < 3; i++){
    pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
    //for(cs_lnum_t j = 0; j < 3; j++){
      fctb[j] = dc[j] * pfac;
      Cell::ref(rhs[c_id1][i][j]).conflict_free_add(-1u, Cell::ref(_weight2 * fctb[j]));
      Cell::ref(rhs[c_id2][i][j]).conflict_free_add(-1u, Cell::ref(_weight1 * fctb[j]));
      //atomicAdd(&rhs[c_id1][i][j], _weight2 * fctb[j]);
      //atomicAdd(&rhs[c_id2][i][j], _weight1 * fctb[j]);
    //}
  //}
}

__global__ static void
_compute_gradient_lsq_v_v5(cs_lnum_t           n_cells,
                        cs_real_t        *restrict _gradv,
                        cs_real_t        *restrict _rhs,
                        cs_cocg_6_t         *restrict cocg)
{
  cs_real_t *rhs = (cs_real_t *) _rhs;
  cs_real_t *gradv = (cs_real_t *) _gradv;
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells) 
    return;

  size_t c_id1 = c_id / (3*3);
  size_t i = (c_id / 3) % 3;
  size_t j = c_id % 3;

  auto cocg_temp = cocg[c_id1];
  cs_real_t _cocg[3];

  _cocg[0] = cocg_temp[5];
  _cocg[1] = cocg_temp[4];
  _cocg[2] = cocg_temp[2];

  if(j == 0){
    _cocg[0] = cocg_temp[0];
    _cocg[1] = cocg_temp[3];
    _cocg[2] = cocg_temp[5];
  }

  if(j == 1){
    _cocg[0] = cocg_temp[3];
    _cocg[1] = cocg_temp[1];
    _cocg[2] = cocg_temp[4];
  }
  
  gradv[c_id] =   rhs[c_id1*3*3 + i*3] * _cocg[0]
                        + rhs[c_id1*3*3 + i*3 + 1] * _cocg[1]
                        + rhs[c_id1*3*3 + i*3 + 2] * _cocg[2];

}

__global__ static void
_compute_gradient_lsq_v_v6(cs_lnum_t           n_cells,
                        cs_real_33_t        *restrict gradv,
                        cs_real_33_t        *restrict rhs,
                        cs_cocg_6_t         *restrict cocg)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells) 
    return;

  size_t c_id1 = c_id / (3*3);
  size_t i = (c_id / 3) % 3;
  size_t j = c_id % 3;

  auto cocg_temp = cocg[c_id1];
  cs_real_t _cocg[3];

  _cocg[0] = cocg_temp[5];
  _cocg[1] = cocg_temp[4];
  _cocg[2] = cocg_temp[2];

  if(j == 0){
    _cocg[0] = cocg_temp[0];
    _cocg[1] = cocg_temp[3];
    _cocg[2] = cocg_temp[5];
  }

  if(j == 1){
    _cocg[0] = cocg_temp[3];
    _cocg[1] = cocg_temp[1];
    _cocg[2] = cocg_temp[4];
  }

  gradv[c_id1][i][j] =   rhs[c_id1][i][0] * _cocg[0]
                        + rhs[c_id1][i][1] * _cocg[1]
                        + rhs[c_id1][i][2] * _cocg[2];

}
