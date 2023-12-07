__global__ static void
cs_slope_test_gradient_vector_cuda_i( const cs_lnum_t             n_i_faces,
                                    const cs_lnum_2_t *restrict   i_face_cells,
                                    const cs_real_3_t *restrict   i_face_cog,
                                    const cs_real_3_t *restrict   cell_cen,
                                    const cs_real_3_t            *pvar,
                                    const cs_real_t *restrict     i_massflux,
                                    const cs_real_3_t *restrict   i_f_face_normal,
                                    cs_real_33_t                 *grad,
                                    cs_real_33_t                 *grdpa)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }
  cs_real_t difv[3], djfv[3], vfac[3];
  cs_real_t pif, pjf, pfac;
  cs_lnum_t c_id1, c_id2;

  c_id1 = i_face_cells[f_id][0];
  c_id2 = i_face_cells[f_id][1];


  for (int jsou = 0; jsou < 3; jsou++) {
    difv[jsou] = i_face_cog[f_id][jsou] - cell_cen[c_id1][jsou];
    djfv[jsou] = i_face_cog[f_id][jsou] - cell_cen[c_id2][jsou];
  }

  /* x-y-z component, p = u, v, w */

  for (int isou = 0; isou < 3; isou++) {
    pif = pvar[c_id1][isou];
    pjf = pvar[c_id2][isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      pif = pif + grad[c_id1][isou][jsou]*difv[jsou];
      pjf = pjf + grad[c_id2][isou][jsou]*djfv[jsou];
    }

    pfac = pjf;
    if (i_massflux[f_id] > 0.) pfac = pif;

    /* U gradient */

    for (int jsou = 0; jsou < 3; jsou++) {
      vfac[jsou] = pfac*i_f_face_normal[f_id][jsou];
      atomicAdd(&grdpa[c_id1][isou][jsou],  vfac[jsou]);
      atomicAdd(&grdpa[c_id2][isou][jsou],- vfac[jsou]);
    }
  }
}


__global__ static void
cs_slope_test_gradient_vector_cuda_b(const cs_lnum_t             n_b_faces,
                                     const cs_real_3_t          *pvar,
                                     const cs_lnum_t *restrict   b_face_cells,
                                     const cs_real_3_t *restrict diipb,
                                     const int                   inc,
                                     const cs_real_3_t          *coefa,
                                     const cs_real_33_t         *coefb,
                                     const cs_real_3_t *restrict b_f_face_normal,
                                     const cs_real_33_t         *grad,
                                     cs_real_33_t               *grdpa)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_b_faces){
    return;
  }

  cs_real_t diipbv[3];
  cs_lnum_t ii = b_face_cells[f_id];

  for (int jsou = 0; jsou < 3; jsou++)
    diipbv[jsou] = diipb[f_id][jsou];

  /* x-y-z components, p = u, v, w */

  for (int isou = 0; isou < 3; isou++) {
    cs_real_t pfac = inc*coefa[f_id][isou];
    /*coefu is a matrix */
    for (int jsou =  0; jsou < 3; jsou++)
      pfac += coefb[f_id][jsou][isou]*(  pvar[ii][jsou]
                                          + grad[ii][jsou][0]*diipbv[0]
                                          + grad[ii][jsou][1]*diipbv[1]
                                          + grad[ii][jsou][2]*diipbv[2]);

    for (int jsou = 0; jsou < 3; jsou++)
      atomicAdd(&grdpa[ii][isou][jsou], pfac*b_f_face_normal[f_id][jsou]);
  }

}



__global__ static void
cs_slope_test_gradient_vector_cuda_f(const cs_lnum_t  n_cells,
                                     cs_real_t       *cell_vol,
                                     cs_real_33_t    *grdpa)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id >= n_cells){
    return;
  }
  size_t c_idx = c_id / (3*3);
  size_t i = (c_id / 3) % 3;
  size_t j = c_id % 3;

  cs_real_t unsvol = 1./cell_vol[c_idx];
  grdpa[c_idx][i][j] *= unsvol;
}