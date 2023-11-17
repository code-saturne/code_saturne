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
_compute_reconstruct_v_i_face_v2(cs_lnum_t            n_i_faces,
                          const cs_lnum_t      *i_group_index,
                          const cs_lnum_2_t      *i_face_cells,
                          const cs_real_3_t    *pvar,
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight,
                          const cs_real_33_t        *restrict r_grad,
                          cs_real_33_t        *restrict grad,
                          const cs_real_3_t *restrict dofij,
                          const cs_real_3_t *restrict i_f_face_normal)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }

  size_t f_idt = f_id / 3;
  size_t i = f_id % 3;

  cs_lnum_t c_id1, c_id2;
  cs_real_t pond, ktpond, pfaci, pfacj, rfac;

  c_id1 = i_face_cells[f_idt][0];
  c_id2 = i_face_cells[f_idt][1];

  pond = weight[f_idt];
  ktpond = (c_weight == NULL) ?
        pond :                    // no cell weighting
        pond * c_weight[c_id1] // cell weighting active
          / (      pond * c_weight[c_id1]
            + (1.0-pond)* c_weight[c_id2]);


  pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
  pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

  /* Reconstruction part */
  rfac = 0.5 * (  dofij[f_idt][0]*(  r_grad[c_id1][i][0]
                                            + r_grad[c_id2][i][0])
                          + dofij[f_idt][1]*(  r_grad[c_id1][i][1]
                                            + r_grad[c_id2][i][1])
                          + dofij[f_idt][2]*(  r_grad[c_id1][i][2]
                                            + r_grad[c_id2][i][2]));

  for (cs_lnum_t j = 0; j < 3; j++) {
    atomicAdd(&grad[c_id1][i][j],(pfaci + rfac) * i_f_face_normal[f_idt][j]);
    atomicAdd(&grad[c_id2][i][j], - ((pfacj + rfac) * i_f_face_normal[f_idt][j]));
  }
    
}




__global__ static void
_compute_reconstruct_v_i_face_v2cf(cs_lnum_t            n_i_faces,
                          const cs_lnum_t      *i_group_index,
                          const cs_lnum_2_t      *i_face_cells,
                          const cs_real_3_t    *pvar,
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight,
                          const cs_real_33_t        *restrict r_grad,
                          cs_real_33_t        *restrict grad,
                          const cs_real_3_t *restrict dofij,
                          const cs_real_3_t *restrict i_f_face_normal)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }

  size_t f_idt = f_id / 3;
  size_t i = f_id % 3;

  cs_lnum_t c_id1, c_id2;
  cs_real_t pond, ktpond, pfaci, pfacj, rfac;

  c_id1 = i_face_cells[f_idt][0];
  c_id2 = i_face_cells[f_idt][1];

  pond = weight[f_idt];
  ktpond = (c_weight == NULL) ?
        pond :                    // no cell weighting
        pond * c_weight[c_id1] // cell weighting active
          / (      pond * c_weight[c_id1]
            + (1.0-pond)* c_weight[c_id2]);


  pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
  pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

  /* Reconstruction part */
  rfac = 0.5 * (  dofij[f_idt][0]*(  r_grad[c_id1][i][0]
                                            + r_grad[c_id2][i][0])
                          + dofij[f_idt][1]*(  r_grad[c_id1][i][1]
                                            + r_grad[c_id2][i][1])
                          + dofij[f_idt][2]*(  r_grad[c_id1][i][2]
                                            + r_grad[c_id2][i][2]));

  using Cell = AtomicCell<cs_real_t,3>;
  Cell grad_cf1, grad_cf2;

  for (cs_lnum_t j = 0; j < 3; j++) {
    grad_cf1[j].get() = (pfaci + rfac) * i_f_face_normal[f_idt][j];
    grad_cf2[j].get() = - ((pfacj + rfac) * i_f_face_normal[f_idt][j]);
    // Cell::ref(grad_cf1[c_id1][i][j]).conflict_free_add(-1u, Cell::ref((pfaci + rfac) * i_f_face_normal[f_idt][j]));
    // Cell::ref(grad_cf2[c_id2][i][j]).conflict_free_add(-1u, Cell::ref(- ((pfacj + rfac) * i_f_face_normal[f_idt][j])));
  }
  Cell::ref(grad[c_id1][i]).conflict_free_add(-1u, grad_cf1);
  Cell::ref(grad[c_id2][i]).conflict_free_add(-1u, grad_cf2);
    
}


__global__ static void
_compute_reconstruct_v_b_face_v2(cs_lnum_t            n_b_faces,
                              const bool                *coupled_faces,
                              cs_lnum_t                 cpl_stride,
                              const cs_real_33_t  *restrict coefbv,
                              const cs_real_3_t   *restrict coefav,
                              const cs_real_3_t   *restrict pvar,
                              int                           inc,
                              const cs_real_3_t *restrict diipb,
                              const cs_real_33_t        *restrict r_grad,
                              cs_real_33_t        *restrict grad,
                              const cs_real_3_t *restrict b_f_face_normal,
                              const cs_lnum_t *restrict b_face_cells)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_b_faces){
    return;
  }

  size_t f_idt = f_id / 3;
  size_t i = f_id % 3;

  cs_lnum_t c_id;
  cs_real_t pond, ktpond, pfac, rfac, vecfac;

  // if (coupled_faces[f_idt * cpl_stride])
  //   return;

  c_id = b_face_cells[f_idt];

  pfac = inc*coefav[f_idt][i];

  for (cs_lnum_t k = 0; k < 3; k++){
    pfac += coefbv[f_idt][i][k] * pvar[c_id][k];
  }

  pfac -= pvar[c_id][i];

//   /* Reconstruction part */
  rfac = 0.;
  for (cs_lnum_t k = 0; k < 3; k++) {
    vecfac =   r_grad[c_id][k][0] * diipb[f_idt][0]
                        + r_grad[c_id][k][1] * diipb[f_idt][1]
                        + r_grad[c_id][k][2] * diipb[f_idt][2];
    rfac += coefbv[f_idt][i][k] * vecfac;
  }

  for (cs_lnum_t j = 0; j < 3; j++){
    atomicAdd(&grad[c_id][i][j], (pfac + rfac) * b_f_face_normal[f_idt][j]);
  }

}



__global__ static void
_compute_reconstruct_v_b_face_v2cf(cs_lnum_t            n_b_faces,
                              const bool                *coupled_faces,
                              cs_lnum_t                 cpl_stride,
                              const cs_real_33_t  *restrict coefbv,
                              const cs_real_3_t   *restrict coefav,
                              const cs_real_3_t   *restrict pvar,
                              int                           inc,
                              const cs_real_3_t *restrict diipb,
                              const cs_real_33_t        *restrict r_grad,
                              cs_real_33_t        *restrict grad,
                              const cs_real_3_t *restrict b_f_face_normal,
                              const cs_lnum_t *restrict b_face_cells)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_b_faces){
    return;
  }

  size_t f_idt = f_id / 3;
  size_t i = f_id % 3;

  cs_lnum_t c_id;
  cs_real_t pond, ktpond, pfac, rfac, vecfac;

  // if (coupled_faces[f_idt * cpl_stride])
  //   return;

  c_id = b_face_cells[f_idt];

  pfac = inc*coefav[f_idt][i];

  for (cs_lnum_t k = 0; k < 3; k++){
    pfac += coefbv[f_idt][i][k] * pvar[c_id][k];
  }

  pfac -= pvar[c_id][i];

//   /* Reconstruction part */
  rfac = 0.;
  for (cs_lnum_t k = 0; k < 3; k++) {
    vecfac =   r_grad[c_id][k][0] * diipb[f_idt][0]
                        + r_grad[c_id][k][1] * diipb[f_idt][1]
                        + r_grad[c_id][k][2] * diipb[f_idt][2];
    rfac += coefbv[f_idt][i][k] * vecfac;
  }

  using Cell = AtomicCell<cs_real_t,3>;
  Cell grad_cf;

  for (cs_lnum_t j = 0; j < 3; j++){
    grad_cf[j].get() = (pfac + rfac) * b_f_face_normal[f_idt][j];
    // grad[c_id][i][j].get() += (pfac + rfac) * b_f_face_normal[f_idt][j];
  }
  Cell::ref(grad[c_id][i]).conflict_free_add(-1u, grad_cf);

}


__global__ static void
_compute_reconstruct_correction_v2(  cs_lnum_t                       n_cells,
                                    cs_lnum_t                       has_dc,
                                    const int *restrict             c_disable_flag,
                                    const cs_real_t *restrict       cell_f_vol,
                                    cs_real_33_t        *restrict   grad,
                                    const cs_real_33_t *restrict    corr_grad_lin,
                                    bool                            test_bool
                                  )
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  

  if(c_id >= n_cells){
    return;
  }
  
  size_t c_idt = c_id / 3;
  size_t i = c_id % 3;

  cs_real_t dvol;
  /* Is the cell disabled (for solid or porous)? Not the case if coupled */
  if (has_dc * c_disable_flag[has_dc * c_idt] == 0)
    dvol = 1. / cell_f_vol[c_idt];
  else
    dvol = 0.;

  for (cs_lnum_t j = 0; j < 3; j++){
    grad[c_idt][i][j] *= dvol;
  }


  if (test_bool) {
    cs_real_t gradpa[3];
    for (cs_lnum_t j = 0; j < 3; j++) {
      gradpa[j] = grad[c_idt][i][j];
    }

    for (cs_lnum_t j = 0; j < 3; j++) {
      grad[c_idt][i][j] = corr_grad_lin[c_idt][j][0] * gradpa[0]
                         + corr_grad_lin[c_idt][j][1] * gradpa[1]
                         + corr_grad_lin[c_idt][j][2] * gradpa[2];
    }
  }

}