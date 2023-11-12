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
_compute_rhs_lsq_v_i_face_gather(cs_lnum_t            n_cells,
                          const cs_lnum_t      *restrict cell_cells_idx,
                          const cs_lnum_t      *restrict cell_cells,
                          const cs_lnum_t      *restrict cell_i_faces,
                          const short int      *restrict cell_i_faces_sgn,
                          const cs_real_3_t    *restrict cell_f_cen,
                          cs_real_33_t         *restrict rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t         *restrict weight,
                          const cs_real_t      *restrict c_weight)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id1 >= n_cells){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _denom, _pond, pfac;
  cs_lnum_t c_id2, f_id;

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];

  for(cs_lnum_t index = s_id; index < e_id; index++){
    c_id2 = cell_cells[index];

    dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
    dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
    dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

    ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
    
    if (c_weight != NULL){
      f_id = cell_i_faces[index];
    _pond = (cell_i_faces_sgn[index] > 0) ? weight[f_id] : 1. - weight[f_id];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
                            
    for(cs_lnum_t i = 0; i < 3; i++){
      pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
      for(cs_lnum_t j = 0; j < 3; j++){
        fctb[j] = dc[j] * pfac;
        rhs[c_id1][i][j] += c_weight[c_id2] * _denom * fctb[j];
      }
    }
  }
  else{
    for(cs_lnum_t i = 0; i < 3; i++){
      pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
      for(cs_lnum_t j = 0; j < 3; j++){
        fctb[j] = dc[j] * pfac;
        rhs[c_id1][i][j] += fctb[j];
      }
    }
  }
}
}

__global__ static void
_compute_rhs_lsq_v_b_face_gather(cs_lnum_t           n_b_cells,
                          const cs_lnum_t      *restrict cell_b_faces_idx,
                          const cs_lnum_t      *restrict cell_b_faces,
                          const cs_lnum_t      *restrict b_cells,
                          const cs_real_3_t    *restrict b_face_normal,
                          cs_real_33_t         *restrict rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t      *restrict b_dist,
                          const cs_real_33_t   *restrict coefbv,
                          const cs_real_3_t    *restrict coefav,
                          const int            inc)
{
  cs_lnum_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_idx >= n_b_cells){
    return;
  }

  cs_lnum_t c_id = b_cells[c_idx];

  cs_lnum_t f_id;
  cs_real_t n_d_dist[3], d_b_dist, pfac, norm, inverse_norm;

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  for(cs_lnum_t index = s_id; index < e_id; index++){

    f_id = cell_b_faces[index];

    cs_math_3_normalise_cuda(b_face_normal[f_id], n_d_dist);

    d_b_dist = 1. / b_dist[f_id];

    /* Normal divided by b_dist */
    n_d_dist[0] *= d_b_dist;
    n_d_dist[1] *= d_b_dist;
    n_d_dist[2] *= d_b_dist;

    for (cs_lnum_t i = 0; i < 3; i++) {
      pfac =   coefav[f_id][i]*inc
            + ( coefbv[f_id][0][i] * pvar[c_id][0]
              + coefbv[f_id][1][i] * pvar[c_id][1]
              + coefbv[f_id][2][i] * pvar[c_id][2]
              - pvar[c_id][i]);

      rhs[c_id][i][0] += n_d_dist[0] * pfac;
      rhs[c_id][i][1] += n_d_dist[1] * pfac;
      rhs[c_id][i][2] += n_d_dist[2] * pfac; 
    }
  }
}
