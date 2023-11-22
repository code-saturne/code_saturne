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

template <cs_lnum_t stride, typename val_t, typename coefb_t>
__global__ static void
_compute_rhs_lsq_v_b_face_gather_stride(cs_lnum_t           n_b_cells,
                          const cs_lnum_t      *restrict cell_b_faces_idx,
                          const cs_lnum_t      *restrict cell_b_faces,
                          const cs_lnum_t      *restrict b_cells,
                          const cs_real_3_t    *restrict b_face_cog,
                          const cs_real_3_t    *restrict cell_cen,
                          cs_real_33_t         *restrict rhs,
                          const val_t    *restrict pvar,
                          const coefb_t   *restrict coefbv,
                          const cs_real_3_t    *restrict coefav,
                          cs_cocg_6_t         *restrict cocg,
                          const cs_cocg_6_t         *restrict cocgb,
                          const int            inc)
{
  cs_lnum_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_idx >= n_b_cells){
    return;
  }

  cs_lnum_t c_id = b_cells[c_idx];

  cs_lnum_t f_id;
  cs_real_t dif[stride], ddif, pfac, norm, var_f[stride];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  for(cs_lnum_t ll = 0; ll < 6; ll++)
    cocg[c_id][ll] = cocgb[c_idx][ll];

  for(cs_lnum_t index = s_id; index < e_id; index++){

    f_id = cell_b_faces[index];

    for (cs_lnum_t ll = 0; ll < 3; ll++)
      dif[ll] = b_face_cog[f_id][ll] - cell_cen[c_id][ll];

    ddif = 1. / cs_math_3_square_norm_cuda(dif);

    cocg[c_id][0] += dif[0]*dif[0]*ddif;
    cocg[c_id][1] += dif[1]*dif[1]*ddif;
    cocg[c_id][2] += dif[2]*dif[2]*ddif;
    cocg[c_id][3] += dif[0]*dif[1]*ddif;
    cocg[c_id][4] += dif[1]*dif[2]*ddif;
    cocg[c_id][5] += dif[0]*dif[2]*ddif;

    for (cs_lnum_t kk = 0; kk < stride; kk++) {
      var_f[kk] = coefav[f_id][kk]*inc;
      for (cs_lnum_t ll = 0; ll < stride; ll++) {
        var_f[kk] += coefbv[f_id][ll][kk] * pvar[c_id][ll];
      }

      pfac = (var_f[kk] - pvar[c_id][kk]) * ddif;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhs[c_id][kk][ll] += dif[ll] * pfac;
    }
  }
  _math_6_inv_cramer_sym_in_place_cuda(cocg[c_id]);
}
