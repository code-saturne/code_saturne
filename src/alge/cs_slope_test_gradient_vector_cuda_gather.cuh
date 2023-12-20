__global__ static void
cs_slope_test_gradient_vector_cuda_i_gather( const cs_lnum_t      n_cells,
                                    const cs_real_3_t *restrict   i_face_cog,
                                    const cs_real_3_t *restrict   cell_cen,
                                    const cs_real_3_t            *pvar,
                                    const cs_real_t *restrict     i_massflux,
                                    const cs_real_3_t *restrict   i_f_face_normal,
                                    const cs_lnum_t *restrict     cell_cells_idx,
                                    const cs_lnum_t *restrict     cell_cells,
                                    const cs_lnum_t *restrict     cell_i_faces,
                                    const short int *restrict     cell_i_faces_sgn,
                                    cs_real_33_t                 *grad,
                                    cs_real_33_t                 *grdpa)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id1 >= n_cells){
    return;
  }

  cs_real_t difv[3], djfv[3], vfac[3];
  cs_real_t pif, pjf, pfac, face_sgn;
  cs_lnum_t c_id2, f_id;

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];


  for(cs_lnum_t index = s_id; index < e_id; index++){
    c_id2 = cell_cells[index];
    f_id = cell_i_faces[index];
    face_sgn = cell_i_faces_sgn[index];

    for (int jsou = 0; jsou < 3; jsou++) {
      difv[jsou] = i_face_cog[f_id][jsou] - cell_cen[c_id1][jsou];
      djfv[jsou] = i_face_cog[f_id][jsou] - cell_cen[c_id2][jsou];
    }

    for (int isou = 0; isou < 3; isou++) {
      pif = pvar[c_id1][isou];
      pjf = pvar[c_id2][isou];
      for (int jsou = 0; jsou < 3; jsou++) {
        pif = pif + grad[c_id1][isou][jsou] * difv[jsou];
        pjf = pjf + grad[c_id2][isou][jsou] * djfv[jsou];
      }

      pfac = pjf;
      if (i_massflux[f_id] * face_sgn > 0.) 
        pfac = pif;
      pfac *= face_sgn;

      for (int jsou = 0; jsou < 3; jsou++) {
        vfac[jsou] = pfac*i_f_face_normal[f_id][jsou];
        grdpa[c_id1][isou][jsou] += vfac[jsou];
      }
    }
  }
}


__global__ static void
cs_slope_test_gradient_vector_cuda_b_gather(const cs_lnum_t      n_b_cells,
                                     const cs_real_3_t          *pvar,
                                     const cs_real_3_t *restrict diipb,
                                     const int                   inc,
                                     const cs_real_3_t          *coefa,
                                     const cs_real_33_t         *coefb,
                                     const cs_real_3_t *restrict b_f_face_normal,
                                     const cs_lnum_t   *restrict b_cells,
                                     const cs_lnum_t   *restrict cell_b_faces,
                                     const cs_lnum_t   *restrict cell_b_faces_idx,
                                     const cs_real_33_t         *grad,
                                     cs_real_33_t               *grdpa)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;


  if(c_id1 >= n_b_cells){
    return;
  }

  cs_lnum_t c_id = b_cells[c_id1];
  
  cs_real_t pfac, rfac, vecfac;
  cs_real_t diipbv[3];
  cs_lnum_t f_id;
  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  for(cs_lnum_t index = s_id; index < e_id; index++){
    f_id = cell_b_faces[index];

    for (int jsou = 0; jsou < 3; jsou++)
      diipbv[jsou] = diipb[f_id][jsou];

    for (int isou = 0; isou < 3; isou++) {
      pfac = inc*coefa[f_id][isou];
      /*coefu is a matrix */
      for (int jsou =  0; jsou < 3; jsou++)
        pfac += coefb[f_id][jsou][isou]*(  pvar[c_id][jsou]
                                            + grad[c_id][jsou][0]*diipbv[0]
                                            + grad[c_id][jsou][1]*diipbv[1]
                                            + grad[c_id][jsou][2]*diipbv[2]);

      for (int jsou = 0; jsou < 3; jsou++)
        grdpa[c_id][isou][jsou] += pfac*b_f_face_normal[f_id][jsou];
    }
  }
}
