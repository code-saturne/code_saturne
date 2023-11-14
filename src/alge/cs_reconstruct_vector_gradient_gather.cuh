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
_compute_reconstruct_v_i_face_gather(cs_lnum_t            n_cells,
                          const cs_lnum_2_t      *i_face_cells,
                          const cs_real_3_t    *pvar,
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight,
                          const cs_real_33_t        *restrict r_grad,
                          cs_real_33_t        *restrict grad,
                          const cs_real_3_t *restrict dofij,
                          const cs_real_3_t *restrict i_f_face_normal,
                          const cs_lnum_t *restrict cell_cells_idx,
                          const cs_lnum_t *restrict cell_cells,
                          const cs_lnum_t *restrict cell_i_faces,
                          const short int *restrict cell_i_faces_sgn,
                          const cs_lnum_t n_i_faces)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id1 >= n_cells){
    return;
  }

  cs_lnum_t c_id2, f_id;
  cs_real_t pond, ktpond, pfaci, pfacj, rfac;

  // if(cell_cells_idx) printf("erreur dans le kernel");
  // if(cell_cells) printf("erreur dans le kernel");
  // if(cell_i_faces) printf("erreur dans le kernel");
  // if(cell_i_faces_sgn) printf("erreur dans le kernel");

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];
  // printf("s_id = %d\t",s_id);
  // printf("e_id = %d\t",e_id);


  for(cs_lnum_t index = s_id; index < e_id; index++){
    c_id2 = cell_cells[index];
    f_id = cell_i_faces[index];

    // pond = weight[f_id];
    pond = (cell_i_faces_sgn[index] > 0) ? weight[f_id] : 1. - weight[f_id];
    ktpond = (c_weight == NULL) ?
            pond :                    // no cell weighting
            pond * c_weight[c_id1] // cell weighting active
            / (      pond * c_weight[c_id1]
                + (1.0-pond)* c_weight[c_id2]);


    for (cs_lnum_t i = 0; i < 3; i++) {
      pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
      // pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

      /* Reconstruction part */
      rfac = 0.5 * (  dofij[f_id][0]*(  r_grad[c_id1][i][0]
                                                + r_grad[c_id2][i][0])
                              + dofij[f_id][1]*(  r_grad[c_id1][i][1]
                                                + r_grad[c_id2][i][1])
                              + dofij[f_id][2]*(  r_grad[c_id1][i][2]
                                                + r_grad[c_id2][i][2]));

      for (cs_lnum_t j = 0; j < 3; j++) {
        grad[c_id1][i][j] += (pfaci + rfac) * i_f_face_normal[f_id][j];
        // grad[c_id1][i][j] -= (pfacj + rfac) * i_f_face_normal[f_id][j];
      }
    }
    // grad[0][0][0] = 1. + f_id;
    // grad[c_id2][0][0] = 1. + pond;
  }
}
