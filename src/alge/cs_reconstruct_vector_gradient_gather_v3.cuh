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


template <cs_lnum_t stride>
__global__ static void
_compute_reconstruct_v_i_face_gather_v3(cs_lnum_t            n_cells,
                          const cs_real_t (*restrict pvar)[stride],
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight,
                          const cs_real_t (*restrict r_grad)[stride][3],
                          cs_real_t (*restrict grad)[stride][3],
                          const cs_real_3_t *restrict dofij,
                          const cs_real_3_t *restrict i_f_face_normal,
                          const cs_lnum_t *restrict cell_cells_idx,
                          const cs_lnum_t *restrict cell_cells,
                          const cs_lnum_t *restrict cell_i_faces,
                          const short int *restrict cell_i_faces_sgn)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id1 >= n_cells){
    return;
  }

  cs_lnum_t c_id2, f_id;
  cs_real_t pond, ktpond, pfaci, pfacj, rfac;

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];
  
  auto _grad = grad[c_id1];
  auto _pvar1 = pvar[c_id1];
  auto _r_grad1 = r_grad[c_id1];

  for(cs_lnum_t index = s_id; index < e_id; index++){
    c_id2 = cell_cells[index];
    f_id = cell_i_faces[index];

    auto _pvar2 = pvar[c_id2];
    auto _r_grad2 = r_grad[c_id2];
    auto _dofij =  dofij[f_id];
    auto _i_f_face_normal =  i_f_face_normal[f_id];
    auto _cell_i_faces_sgn =  cell_i_faces_sgn[index];

    pond = (_cell_i_faces_sgn > 0) ? weight[f_id] : 1. - weight[f_id];
    ktpond = (c_weight == NULL) ?
            pond :                    // no cell weighting
            pond * c_weight[c_id1] // cell weighting active
            / (      pond * c_weight[c_id1]
                + (1.0-pond)* c_weight[c_id2]);

    for (cs_lnum_t i = 0; i < stride; i++) {
        pfaci = (1.0-ktpond) * (_pvar2[i] - _pvar1[i]);

        /* Reconstruction part */
        rfac = 0.5 * (    _dofij[0]*(       _r_grad1[i][0]
                                          + _r_grad2[i][0])
                        + _dofij[1]*(       _r_grad1[i][1]
                                          + _r_grad2[i][1])
                        + _dofij[2]*(       _r_grad1[i][2]
                                          + _r_grad2[i][2]));

        for (cs_lnum_t j = 0; j < 3; j++) {
            _grad[i][j] += _cell_i_faces_sgn * (pfaci + rfac) * _i_f_face_normal[j];
        }
    }
  }
  for(cs_lnum_t i = 0; i < stride; i++){
    for(cs_lnum_t j = 0; j < 3; j++){
      grad[c_id1][i][j] = _grad[i][j];
    }
  }
}



template <cs_lnum_t stride>
__global__ static void
_compute_reconstruct_v_b_face_gather_v3(cs_lnum_t           n_b_cells,
                              const cs_real_t (*restrict coefbv)[stride][stride],
                              const cs_real_t (*restrict coefav)[stride],
                              const cs_real_t (*restrict pvar)[stride],
                              int                           inc,
                              const cs_real_3_t *restrict diipb,
                              const cs_real_t (*restrict r_grad)[stride][3],
                              cs_real_t (*restrict grad)[stride][3],
                              const cs_real_3_t *restrict b_f_face_normal,
                              const cs_lnum_t      *restrict b_cells,
                              const cs_lnum_t      *restrict cell_b_faces,
                              const cs_lnum_t      *restrict cell_b_faces_idx)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;


  if(c_id1 >= n_b_cells){
    return;
  }

  cs_lnum_t c_id = b_cells[c_id1];
  
  cs_real_t pfac, rfac, vecfac;
  cs_lnum_t f_id;
  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  auto _grad = grad[c_id];
  auto _r_grad = r_grad[c_id];
  auto _pvar = pvar[c_id];

  for(cs_lnum_t index = s_id; index < e_id; index++){
    f_id = cell_b_faces[index];

    auto _diipb = diipb[f_id];
    auto _coefav = coefav[f_id];
    auto _coefbv = coefbv[f_id];
    auto _b_f_face_normal = b_f_face_normal[f_id];

    for (cs_lnum_t i = 0; i < stride; i++) {

      pfac = inc*_coefav[i];

      for (cs_lnum_t k = 0; k < 3; k++){
        pfac += _coefbv[i][k] * _pvar[k];
      }

      pfac -= _pvar[i];

    //   /* Reconstruction part */
      rfac = 0.;
      for (cs_lnum_t k = 0; k < stride; k++) {
        vecfac =   _r_grad[k][0] * _diipb[0]
                 + _r_grad[k][1] * _diipb[1]
                 + _r_grad[k][2] * _diipb[2];
        rfac += _coefbv[i][k] * vecfac;
      }

      for (cs_lnum_t j = 0; j < 3; j++){
        _grad[i][j] += (pfac + rfac) * _b_f_face_normal[j];
      }

    }
  }
  for(cs_lnum_t i = 0; i < stride; i++){
    for(cs_lnum_t j = 0; j < 3; j++){
      grad[c_id][i][j] = _grad[i][j];
    }
  }
  
}
