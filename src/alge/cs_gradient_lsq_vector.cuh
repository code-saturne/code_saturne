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
_init_rhs(cs_lnum_t         n_cells_ext,
           cs_real_33_t      *restrict rhs)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id < n_cells_ext) {
    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
        rhs[c_id][i][j] = 0.0;
  }
}

__global__ static void
_compute_rhs_lsq_v_i_face_v0(cs_lnum_t            n_i_faces,
                          const cs_lnum_2_t      *i_face_cells,
                          const cs_real_3_t    *cell_f_cen,
                          cs_real_33_t         *rhs,
                          const cs_real_3_t    *pvar,
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= n_i_faces){
    return;
  }
  cs_real_t dc[3], fctb[3], ddc, _weight1, _weight2, _denom, _pond, pfac;
  cs_lnum_t c_id1, c_id2;

  c_id1 = i_face_cells[f_id][0];
  c_id2 = i_face_cells[f_id][1];

  dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
  dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
  dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

  ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);
  
  if (c_weight != NULL){
    _pond = weight[f_id];
    _denom = 1. / (  _pond       *c_weight[c_id1]
                                + (1. - _pond)*c_weight[c_id2]);
                            
    for(cs_lnum_t i = 0; i < 3; i++){
      pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
      for(cs_lnum_t j = 0; j < 3; j++){
        fctb[j] = dc[j] * pfac;
        atomicAdd(&rhs[c_id1][i][j], c_weight[c_id2] * _denom * fctb[j]);
        atomicAdd(&rhs[c_id2][i][j], c_weight[c_id1] * _denom * fctb[j]);
      }
    }
  }
  else{
    for(cs_lnum_t i = 0; i < 3; i++){
      pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;
      for(cs_lnum_t j = 0; j < 3; j++){
        fctb[j] = dc[j] * pfac;
        atomicAdd(&rhs[c_id1][i][j], fctb[j]);
        atomicAdd(&rhs[c_id2][i][j], fctb[j]);
      }
    }
  }
}

__global__ static void
_compute_rhs_lsq_v_i_face(cs_lnum_t            n_i_faces,
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

  c_id1 = i_face_cells[f_id][0];
  c_id2 = i_face_cells[f_id][1];

  dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
  dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
  dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

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
      atomicAdd(&rhs[c_id1][i][j], _weight2 * fctb[j]);
      atomicAdd(&rhs[c_id2][i][j], _weight1 * fctb[j]);
    }
  }
}

__global__ static void
_compute_rhs_lsq_v_i_face_cf(cs_lnum_t            size,
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

  c_id1 = i_face_cells[f_id][0];
  c_id2 = i_face_cells[f_id][1];

  dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
  dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
  dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

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
      _rhs1[i][j].get() = _weight2 * fctb[j];
      _rhs2[i][j].get() = _weight1 * fctb[j];
      //atomicAdd(&rhs[c_id1][i][j], _weight2 * fctb[j]);
      //atomicAdd(&rhs[c_id2][i][j], _weight1 * fctb[j]);
    }
  }

#if 1
  Cell::ref(rhs[c_id1]).conflict_free_add(-1u, _rhs1);
  Cell::ref(rhs[c_id2]).conflict_free_add(-1u, _rhs2);
#else
  Cell::ref(rhs[c_id1]).atomic_add(_rhs1);
  Cell::ref(rhs[c_id2]).atomic_add(_rhs2);
#endif
}

__global__ static void
_compute_rhs_lsq_v_b_neighbor(cs_lnum_t            n_cells,
                                const cs_lnum_t      *restrict cell_cells_idx,
                                const cs_lnum_t      *restrict cell_cells,
                                const cs_real_3_t    *restrict cell_f_cen,
                                cs_real_33_t         *restrict rhs,
                                const cs_real_3_t    *restrict pvar)
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id1 >= n_cells){
    return;
  }

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];

  cs_real_t dc[3], ddc, pfac;

  for(cs_lnum_t index = s_id; index < e_id; index++){

    cs_lnum_t c_id2 = cell_cells[index];

    dc[0] = cell_f_cen[c_id2][0] - cell_f_cen[c_id1][0];
    dc[1] = cell_f_cen[c_id2][1] - cell_f_cen[c_id1][1];
    dc[2] = cell_f_cen[c_id2][2] - cell_f_cen[c_id1][2];

    ddc = 1./(dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]);

    for (cs_lnum_t i = 0; i < 3; i++) {

      pfac = (pvar[c_id2][i] - pvar[c_id1][i]) * ddc;

      for (cs_lnum_t j = 0; j < 3; j++) {
        rhs[c_id1][i][j] += dc[j] * pfac;
      }
    }
  }

}

__global__ static void
_compute_rhs_lsq_v_b_face(cs_lnum_t           n_b_faces,
                          const cs_lnum_t      *restrict b_face_cells,
                          const cs_real_3_t    *restrict cell_f_cen,
                          const cs_real_3_t    *restrict b_face_normal,
                          cs_real_33_t         *restrict rhs,
                          const cs_real_3_t    *restrict pvar,
                          const cs_real_t      *restrict b_dist,
                          const cs_real_33_t   *restrict coefbv,
                          const cs_real_3_t    *restrict coefav,
                          const int            inc)
{
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

    atomicAdd(&rhs[c_id1][i][0], n_d_dist[0] * pfac);
    atomicAdd(&rhs[c_id1][i][1], n_d_dist[1] * pfac);
    atomicAdd(&rhs[c_id1][i][2], n_d_dist[2] * pfac); 
  }
}

__global__ static void
_compute_gradient_lsq_v(cs_lnum_t           n_cells,
                        cs_real_33_t        *restrict gradv,
                        cs_real_33_t        *restrict rhs,
                        cs_cocg_6_t         *restrict cocg)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= n_cells) 
    return;

  for(cs_lnum_t i = 0; i < 3; i++){
    gradv[c_id][i][0] =   rhs[c_id][i][0] * cocg[c_id][0]
                          + rhs[c_id][i][1] * cocg[c_id][3]
                          + rhs[c_id][i][2] * cocg[c_id][5];

    gradv[c_id][i][1] =   rhs[c_id][i][0] * cocg[c_id][3]
                        + rhs[c_id][i][1] * cocg[c_id][1]
                        + rhs[c_id][i][2] * cocg[c_id][4];

    gradv[c_id][i][2] =   rhs[c_id][i][0] * cocg[c_id][5]
                        + rhs[c_id][i][1] * cocg[c_id][4]
                        + rhs[c_id][i][2] * cocg[c_id][2];
  }
}

__global__ static void
_compute_gradient_lsq_b_v(cs_lnum_t           n_b_cells,
                        const cs_lnum_t      *restrict b_cells,
                        const cs_lnum_t      *restrict cell_b_faces_idx,
                        const cs_lnum_t      *restrict cell_b_faces,
                        const cs_real_3_t    *restrict b_face_normal,
                        const cs_real_3_t    *restrict diipb,
                        const cs_real_3_t    *restrict pvar,
                        const cs_real_t      *restrict b_dist,
                        const cs_real_33_t   *restrict coefbv,
                        const cs_real_3_t    *restrict coefav,
                        cs_real_33_t        *restrict gradv,
                        cs_real_33_t        *restrict rhs,
                        cs_cocg_6_t         *restrict cocgb_s,
                        const int            inc)
{
  size_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_idx >= n_b_cells) 
    return;

  cs_lnum_t c_id = b_cells[c_idx];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

  cs_lnum_t f_id;
  cs_real_t cocgb[3][3], cocgb_v[45], rhsb_v[9], x[9];
  cs_real_3_t normal;

  cs_lnum_t _33_9_idx[9][2];
  int nn = 0;
  for (int ll = 0; ll < 3; ll++) {
    for (int mm = 0; mm < 3; mm++) {
      _33_9_idx[nn][0] = ll;
      _33_9_idx[nn][1] = mm;
      nn++;
    }
  }

  auto _cocg = cocgb_s[c_idx];
  auto _rhs = rhs[c_id];

  cocgb[0][0] = _cocg[0];
  cocgb[0][1] = _cocg[3];
  cocgb[0][2] = _cocg[5];
  cocgb[1][0] = _cocg[3];
  cocgb[1][1] = _cocg[1];
  cocgb[1][2] = _cocg[4];
  cocgb[2][0] = _cocg[5];
  cocgb[2][1] = _cocg[4];
  cocgb[2][2] = _cocg[2];

  for(cs_lnum_t index = s_id; index < e_id; index++){
    f_id = cell_b_faces[index];

    cs_math_3_normalize_cuda(b_face_normal[f_id], normal);
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++)
        cocgb[ii][jj] += normal[ii] * normal[jj];
    }
  }

  for (int ll = 0; ll < 9; ll++) {

    int ll_9 = ll*(ll+1)/2;

    for (int mm = 0; mm <= ll; mm++) {
      cocgb_v[ll_9+mm] = 0.;

      int pp = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      int rr = _33_9_idx[mm][0];
      int ss = _33_9_idx[mm][1];

      if (pp == rr)
        cocgb_v[ll_9+mm] = cocgb[qq][ss];

      rhsb_v[ll] = _rhs[pp][qq];
    }
  }

  cs_real_3_t nb;
  cs_real_t a[3], bt[3][3], db, db2;
  for (cs_lnum_t i = s_id; i < e_id; i++) {

    f_id = cell_b_faces[i];

    auto iipbf = diipb[f_id];
    
    cs_math_3_normalize_cuda(b_face_normal[f_id], nb);

    db = 1./b_dist[f_id];
    db2 = db*db;

    for (int ll = 0; ll < 3; ll++) {
      for (int pp = 0; pp < 3; pp++)
        bt[ll][pp] = coefbv[f_id][ll][pp];
    }
    for (int ll = 0; ll < 3; ll++) {
      a[ll] = inc*coefav[f_id][ll];
      bt[ll][ll] -= 1;
    }

    for (int ll = 0; ll < 9; ll++) {

      int kk = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      int ll_9 = ll*(ll+1)/2;
      for (int pp = 0; pp <= ll; pp++) {

        int rr = _33_9_idx[pp][0];
        int ss = _33_9_idx[pp][1];

        cs_real_t cocgv = 0.;
        for (int mm = 0; mm < 3; mm++)
          cocgv += bt[mm][kk]*bt[mm][rr];
        cocgb_v[ll_9+pp] += cocgv*(iipbf[qq]*iipbf[ss])*db2;

        cocgb_v[ll_9+pp] -= (  nb[ss]*bt[rr][kk]*iipbf[qq]
                             + nb[qq]*bt[kk][rr]*iipbf[ss])
                             *db;
      }
    }

    for (int ll = 0; ll < 9; ll++) {
      int pp = _33_9_idx[ll][0];
      int qq = _33_9_idx[ll][1];

      cs_real_t rhsv = 0.;
      for (int rr = 0; rr < 3; rr++) {
        rhsv +=   bt[rr][pp]*diipb[f_id][qq]
                            *(a[rr]+ bt[rr][0]*pvar[c_id][0]
                                   + bt[rr][1]*pvar[c_id][1]
                                   + bt[rr][2]*pvar[c_id][2]);
      }

      rhsb_v[ll] -= rhsv*db2;
    }

  }
  _fact_crout_pp_cuda<9>(cocgb_v);

  _fw_and_bw_ldtl_pp_cuda<9>(cocgb_v, x, rhsb_v);

  for (int kk = 0; kk < 9; kk++) {
    int ii = _33_9_idx[kk][0];
    int jj = _33_9_idx[kk][1];
    gradv[c_id][ii][jj] = x[kk];
  }
}

template <cs_lnum_t stride>
__global__ static void
_compute_gradient_lsq_b_strided_v(const cs_lnum_t           n_b_cells,
                        const cs_lnum_t      *restrict b_cells,
                        const cs_lnum_t      *restrict cell_b_faces_idx,
                        const cs_lnum_t      *restrict cell_b_faces,
                        const cs_real_3_t          *restrict b_face_cog,
                        const cs_real_3_t          *restrict cell_cen,
                        const cs_real_3_t          *restrict diipb,
                        cs_real_t (*restrict gradv)[stride][3],
                        const cs_real_t (*restrict coefbv)[stride][stride],
                        cs_cocg_6_t         *restrict cocg,
                        cs_lnum_t           n_c_iter_max,
                        cs_real_t           c_eps)
{
  size_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_idx >= n_b_cells) 
    return;
  
  cs_lnum_t c_id = b_cells[c_idx];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

  auto c_grad = gradv[c_id];
  cs_real_t grad_0[stride][3], grad_i[stride][3], rhs_c[stride][3], dif[3],  grad_c[stride][3], 
            var_ip_f[stride];

  cs_real_t ref_norm = 0.0, ddif, c_norm = 0;
  cs_lnum_t n_c_it, f_id;
  cs_real_t  eps_dvg = 1e-2;
  cs_real_t cs_math_epzero = 1e-12;

  for(cs_lnum_t i = 0; i < stride; i++){
    for(cs_lnum_t j = 0; j < 3; j++){
      grad_0[i][j] = c_grad[i][j];
      grad_i[i][j] = c_grad[i][j];
    }
  }

  ref_norm = 0;
  for (cs_lnum_t kk = 0; kk < stride; kk++) {
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      ref_norm += cs_math_fabs_cuda(c_grad[kk][ll]);
  }

  c_norm = 0;

  for (n_c_it = 0; n_c_it < n_c_iter_max; n_c_it++) {

    for (cs_lnum_t ll = 0; ll < stride; ll++) {
      rhs_c[ll][0] = 0;
      rhs_c[ll][1] = 0;
      rhs_c[ll][2] = 0;
    }
    
    for(cs_lnum_t index = s_id; index < e_id; index++){
      f_id = cell_b_faces[index];

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        dif[ii] = b_face_cog[f_id][ii] - cell_cen[c_id][ii];

      ddif = 1. / cs_math_3_square_norm_cuda(dif);

      for (cs_lnum_t ll = 0; ll < stride; ll++) {
        var_ip_f[ll] = cs_math_3_dot_product_cuda(c_grad[ll], diipb[f_id]);
      }        

      auto b = coefbv[f_id];

      for (cs_lnum_t kk = 0; kk < stride; kk++) {
        cs_real_t pfac = 0;
        for (cs_lnum_t ll = 0; ll < stride; ll++) {
          pfac += b[kk][ll] * var_ip_f[ll] * ddif;
        }

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhs_c[kk][ll] += dif[ll] * pfac;
      }

    }

    for(cs_lnum_t i = 0; i < stride; i++){
      grad_c[i][0] =  rhs_c[i][0] * cocg[c_id][0]
                    + rhs_c[i][1] * cocg[c_id][3]
                    + rhs_c[i][2] * cocg[c_id][5];

      grad_c[i][1] =  rhs_c[i][0] * cocg[c_id][3]
                    + rhs_c[i][1] * cocg[c_id][1]
                    + rhs_c[i][2] * cocg[c_id][4];

      grad_c[i][2] =  rhs_c[i][0] * cocg[c_id][5]
                    + rhs_c[i][1] * cocg[c_id][4]
                    + rhs_c[i][2] * cocg[c_id][2];
    }

    c_norm = 0.0;
    for (cs_lnum_t ii = 0; ii < stride; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        c_grad[ii][jj] = grad_0[ii][jj] + grad_c[ii][jj];
        c_norm += cs_math_fabs_cuda(c_grad[ii][jj] - grad_i[ii][jj]);
        grad_i[ii][jj] = c_grad[ii][jj];
      }
    }

    if (c_norm < ref_norm * c_eps || c_norm < cs_math_epzero)
        break;
  }
  
  for (cs_lnum_t ii = 0; ii < stride; ii++) {
    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      gradv[c_id][ii][jj] = c_grad[ii][jj];
    }
  }

  if (c_norm > eps_dvg * ref_norm) {
    for (cs_lnum_t ii = 0; ii < stride; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        gradv[c_id][ii][jj] = grad_0[ii][jj];
      }
    }

    n_c_it *= -1;
  }
}
