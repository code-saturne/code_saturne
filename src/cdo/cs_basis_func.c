/*============================================================================
 * Build a set of basis functions for cells and faces and cell gradients
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_log.h"
#include "cs_quadrature.h"
#include "cs_scheme_geometry.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_basis_func.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_BASIS_FUNC_DBG 0

#define CK1_SIZE  4
#define CK2_SIZE  10
#define CKA_SIZE  84  /* Pre-allocated up to order 6 */
#define FK1_SIZE  3
#define FK2_SIZE  6
#define FKA_SIZE  28  /* Pre-allocated up to order 6 */

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double  _clean_threshold = 1e-15;

/* Handle options for building basis functions when using HHO schemes */
static  cs_flag_t  cs_basis_func_hho_face_flag = 0;
static  cs_flag_t  cs_basis_func_hho_cell_flag = 0;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the x^exp where exp is a positive integer
 *
 * \param[in]     x        base
 * \param[in]     exp      exponent
 *
 * \return the result
 */
/*----------------------------------------------------------------------------*/

static inline void
_swap_shint(short int  *a,
            short int  *b)
{
  const short int  tmp = b[0];
  b[0] = a[0], a[0] = tmp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the x^exp where exp is a positive integer
 *
 * \param[in]     x        base
 * \param[in]     exp      exponent
 *
 * \return the result
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_int_pow(cs_real_t    x,
         short int    exp)
{
  assert(exp > -1);
  switch (exp) {

  case 0: /* Unity */
    return 1;
    break;
  case 1: /* Linear */
    return x;
    break;
  case 2: /* Quadratic */
    return x*x;
    break;
  case 3: /* Cubic */
    return x*x*x;
    break;

  default:
    {
      cs_real_t ret = 1.;
      while (exp > 1) {
        if (exp%2 == 0)
          x *= x, exp /= 2;
        else
          ret *= x, x *= x, exp = (exp-1) / 2;
      }
      return x*ret;
    }

  } /* End of switch */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if a 3x3 tensor is nearly diagonal. Nearly means up to a
 *         tolerance given as a parameter.
 *
 * \param[in]  M          3x3 tensor to test
 * \param[in]  tolerance  threshold to consider the tensor as diagonal
 *
 * \return true if tensor is nearly diagonal otherwise false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_is_nearly_diag(const cs_real_t    M[3][3],
                const cs_real_t    tolerance)
{
  assert(fabs(M[0][0]) > cs_math_zero_threshold &&
         fabs(M[1][1]) > cs_math_zero_threshold &&
         fabs(M[2][2]) > cs_math_zero_threshold);

  if (fabs(M[0][1]) + fabs(M[0][2]) > tolerance * fabs(M[0][0]))
    return false;
  if (fabs(M[1][0]) + fabs(M[1][2]) > tolerance * fabs(M[1][1]))
    return false;
  if (fabs(M[2][1]) + fabs(M[2][0]) > tolerance * fabs(M[2][2]))
    return false;

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specific matrix-matrix product R*Q where R is stored in compact way
 *         and Q is defined by its transposed matrix
 *         Case of a 3x3 matrix
 *
 * \param[in]      Qt     transposed of matrix Q
 * \param[in]      R      vector of the coefficient of the decomposition
 * \param[in, out] m      matrix values to compute
 *
 * \Note: R is an upper triangular matrix. Stored in a compact way.
 *
 *  j=   0, 1, 2
 *  i=0| 0| 1| 2|
 *  i=1   | 4| 5|
 *  i=2        6|
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
_product_rq(const cs_real_t    Qt[9],
            const cs_real_t    R[6],
            cs_real_t          m[9])
{
  /* Sanity checks */
  assert(m != NULL && Qt != NULL && R != NULL);

  m[0] = R[0]*Qt[0] + R[1]*Qt[1] + R[2]*Qt[2];  /* M[0][0] */
  m[1] = R[0]*Qt[3] + R[1]*Qt[4] + R[2]*Qt[5];  /* M[0][1] */
  m[2] = R[0]*Qt[6] + R[1]*Qt[7] + R[2]*Qt[8];  /* M[0][2] */
  m[3] =              R[3]*Qt[1] + R[4]*Qt[2];  /* M[1][0] */
  m[4] =              R[3]*Qt[4] + R[4]*Qt[5];  /* M[1][1] */
  m[5] =              R[3]*Qt[7] + R[4]*Qt[8];  /* M[1][2] */
  m[6] =                           R[5]*Qt[2];  /* M[2][0] */
  m[7] =                           R[5]*Qt[5];  /* M[2][1] */
  m[8] =                           R[5]*Qt[8];  /* M[2][2] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the product of Q
 *         and Q is defined by its transposed matrix
 *         Case of a 3x3 matrix.
 *
 * \param[in, out] prodq  matrix values to compute
 * \param[in]      Qt     transposed of matrix Q
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_prodq(cs_real_t          prodq[9],
              const cs_real_t    Qt[9])
{
  /* Sanity checks */
  assert(prodq != NULL && Qt != NULL);
  cs_real_3_t  row = {prodq[0], prodq[1], prodq[2]};

  prodq[0] = row[0]*Qt[0] + row[1]*Qt[1] + row[2]*Qt[2];
  prodq[1] = row[0]*Qt[3] + row[1]*Qt[4] + row[2]*Qt[5];
  prodq[2] = row[0]*Qt[6] + row[1]*Qt[7] + row[2]*Qt[8];

  row[0] = prodq[3], row[1] = prodq[4], row[2] = prodq[5];
  prodq[3] = row[0]*Qt[0] + row[1]*Qt[1] + row[2]*Qt[2];
  prodq[4] = row[0]*Qt[3] + row[1]*Qt[4] + row[2]*Qt[5];
  prodq[5] = row[0]*Qt[6] + row[1]*Qt[7] + row[2]*Qt[8];

  row[0] = prodq[6], row[1] = prodq[7], row[2] = prodq[8];
  prodq[6] = row[0]*Qt[0] + row[1]*Qt[1] + row[2]*Qt[2];
  prodq[7] = row[0]*Qt[3] + row[1]*Qt[4] + row[2]*Qt[5];
  prodq[8] = row[0]*Qt[6] + row[1]*Qt[7] + row[2]*Qt[8];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project a 3-dimension array on another 3-dimension array
 *
 * \param[in, out] a_proj   projection of a_proj on a_ref
 * \param[in]      a_ref    reference
 */
/*----------------------------------------------------------------------------*/

static inline void
_proj3(cs_real_3_t           a_proj,
       const cs_real_3_t     a_ref)
{
  cs_real_t  dp = _dp3(a_proj, a_ref);
  for (int k = 0; k < 3; k++)
    a_proj[k] -= dp * a_ref[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a contribution stemming from a tetrahedron/triangle to the
 *         projection matrix for each quadrature point
 *
 * \param[in]      n_gpts     number of Gauss integration points
 * \param[in]      gpts       list of Gauss points
 * \param[in]      gw         list of weights related to each Gauss point
 * \param[in, out] bf         pointer to a cs_basis_func_t structure
 * \param[in]      first_row  start filling the values from first row
 * \param[in]      n_rows     number of basis functions (n_rows in the matrix)
 * \param[in, out] phi_eval   array storing the evaluation of basis functions
 * \param[in, out] values     values related to the projection matrix
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_contrib(const int           n_gpts,
             const cs_real_3_t   gpts[],
             const cs_real_t     gw[],
             cs_basis_func_t    *bf,
             const int           first_row,
             const int           n_rows,
             cs_real_t           phi_eval[],
             cs_real_t           values[])
{
  for (short int p = 0; p < n_gpts; p++) { /* Loop on Gauss points */

    bf->eval_all_at_point(bf, gpts[p], phi_eval);
    const cs_real_t  w = gw[p];

    /* First row is already computed if k=1 */
    for (short int i = first_row; i < n_rows; i++) {

      const cs_real_t  wphi_i = w * phi_eval[i];
      if (fabs(wphi_i) > cs_math_zero_threshold) {
        cs_real_t  *val_i = values + n_rows*i;
        for (int j = i; j < n_rows; j++)
          val_i[j] += wphi_i * phi_eval[j];
      } /* wphi_i not zero */

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set remaining coefficients of a matrix to get a symmetric one.
 *         Set to zero values below a given threshold
 *
 * \param[in]      n_rows    numer of basis functions (n_rows in the matrix)
 * \param[in, out] values    values related to the projection matrix
 */
/*----------------------------------------------------------------------------*/

static inline void
_symmetrize_and_clean(const int           n_rows,
                      cs_real_t           values[])
{
  for (short int i = 0; i < n_rows; i++) {
    const cs_real_t  *const val_i = values + n_rows*i;
    assert(val_i[i] > 0);
    const cs_real_t  coef = 1/val_i[i];
    for (short int j = i + 1; j < n_rows; j++) {
      if (fabs(val_i[j]*coef) > _clean_threshold)
        values[n_rows*j + i] = val_i[j];
      else
        values[n_rows*j+i] = values[i*n_rows + j] = 0.;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the cell-basis function accordingly to the requested type
 *
 * \param[in, out]  pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]       cm      pointer to a cs_cell_mesh_t
 * \param[in]       id      id of the element to consider
 * \param[in]       center  point used for centering the set of basis functions
 * \param[in, out]  cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_mono_cell_basis_setup(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id,
                       const cs_real_t          center[3],
                       cs_cell_builder_t       *cb)
{
  CS_UNUSED(id);
  CS_UNUSED(cb);

  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  const cs_real_t  odiam = 1. / cm->diam_c;

  for (int k = 0; k < 3; k++) {
    bf->axis[k].meas = odiam;
    bf->center[k] = center[k];
  }

  bf->axis[0].unitv[0] = bf->axis[1].unitv[1] = bf->axis[2].unitv[2] = 1;
  bf->axis[0].unitv[1] = bf->axis[1].unitv[0] = bf->axis[2].unitv[0] = 0;
  bf->axis[0].unitv[2] = bf->axis[1].unitv[2] = bf->axis[2].unitv[1] = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the cell-basis function accordingly to the requested type
 *         Compute the inertia matrix related to the current cell using a QR
 *         algorithm based on a modified Gram-Schmidt algorithm. This is not
 *         an optimal choice but this the matrix is very small (3x3) there is
 *         no real impact
 *
 * \param[in, out]  pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]       cm      pointer to a cs_cell_mesh_t
 * \param[in]       id      id of the element to consider
 * \param[in]       center  point used for centering the set of basis functions
 * \param[in, out]  cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_iner_cell_basis_setup(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id,
                       const cs_real_t          center[3],
                       cs_cell_builder_t       *cb)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ |
                      CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ));

  /* Advanced parameters for controlling the algorithm */
  const int  n_algo_iters = 20;
  const double  tolerance = 1e-12;

  /* Default initialization */
  _mono_cell_basis_setup(bf, cm, id, center, cb);

  cs_real_33_t  M;
  cs_compute_inertia_tensor(cm, center, M);

  if (_is_nearly_diag((const cs_real_t(*)[3])M, tolerance))
    return;  /* Current set of axis are relevant */

  /* QR iterative algorithm to compute eigenvalues and eigenvectors */
  cs_real_6_t  R_k = {0, 0, 0, 0, 0, 0};
  cs_real_9_t  prodQ = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  cs_real_9_t  Qt_k = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  cs_real_9_t  M_k = { M[0][0], M[0][1], M[0][2],
                       M[1][0], M[1][1], M[1][2],
                       M[2][0], M[2][1], M[2][2] };

  int iter;
  for (iter = 0;
       iter < n_algo_iters && !_is_nearly_diag((const cs_real_t(*)[3])M_k,
                                               tolerance);
       iter++) {

    /* Perform the QR factorization using a Modified Gram-Schmidt algo. */
    cs_sdm_33_sym_qr_compute(M_k, Qt_k, R_k);

    /* Compute the product of the successive: prodQ_(k+1) = prodQ_k*Qt_k */
    _update_prodq(prodQ, Qt_k);

    /* M_(k+1) = R_k*Q_k */
    _product_rq(Qt_k, R_k, M_k);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_BASIS_FUNC_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s : cell_id: %d; n_iters: %d\n",
                __func__, cm->c_id, iter);
  cs_log_printf(CS_LOG_DEFAULT,
                "    |% 5.3e % 5.3e % 5.3e|       |% 5.3e % 5.3e % 5.3e|\n"
                " Mk |% 5.3e % 5.3e % 5.3e| ProdQ |% 5.3e % 5.3e % 5.3e|\n"
                "    |% 5.3e % 5.3e % 5.3e|       |% 5.3e % 5.3e % 5.3e|\n",
                M_k[0], M_k[1], M_k[2], prodQ[0], prodQ[1], prodQ[2],
                M_k[3], M_k[4], M_k[5], prodQ[3], prodQ[4], prodQ[5],
                M_k[6], M_k[7], M_k[8], prodQ[6], prodQ[7], prodQ[8]);
#endif

  /* Bubble sort : looking for the greatest eigenvalue */
  const cs_real_3_t  lambdas = {M_k[0], M_k[4], M_k[8]};
  short int  idx[3] = {0, 1, 2};
  if(lambdas[0] < lambdas[1])           _swap_shint(idx,     idx + 1);
  if(lambdas[idx[0]] < lambdas[idx[2]]) _swap_shint(idx,     idx + 2);
  if(lambdas[idx[1]] < lambdas[idx[2]]) _swap_shint(idx + 1, idx + 2);

  /* The eigenvectors are the columns of prodQ */
  cs_real_3_t  a[3];
  a[0][0] = prodQ[idx[0]], a[0][1] = prodQ[idx[0]+3], a[0][2] = prodQ[idx[0]+6];
  a[1][0] = prodQ[idx[1]], a[1][1] = prodQ[idx[1]+3], a[1][2] = prodQ[idx[1]+6];
  a[2][0] = prodQ[idx[2]], a[2][1] = prodQ[idx[2]+3], a[2][2] = prodQ[idx[2]+6];

  /* One final MGS to insure orthogonalities */
  _proj3(a[1], a[0]);
  _proj3(a[2], a[0]);
  _proj3(a[2], a[1]);

  for (int k = 0; k < 3; k++)
    cs_nvec3(a[k], bf->axis + k);

  const cs_real_t  odiam = 1. / cm->diam_c;
  for (int k = 0; k < 3; k++)
    bf->axis[k].meas *= odiam;

#if defined(DEBUG) && !defined(NDEBUG) && CS_BASIS_FUNC_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s : cell_id: %d\n", __func__, cm->c_id);
  cs_log_printf(CS_LOG_DEFAULT,
                "                   unitv                      meas\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n"
                " cell axis |% 5.3e % 5.3e % 5.3e|  % 5.3e\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n",
                bf->axis[0].unitv[0], bf->axis[0].unitv[1],
                bf->axis[0].unitv[2], bf->axis[0].meas,
                bf->axis[1].unitv[0], bf->axis[1].unitv[1],
                bf->axis[1].unitv[2], bf->axis[1].meas,
                bf->axis[2].unitv[0], bf->axis[2].unitv[1],
                bf->axis[2].unitv[2], bf->axis[2].meas);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the face-basis function accordingly to the requested type
 *
 * \param[in, out]  pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]       cm      pointer to a cs_cell_mesh_t
 * \param[in]       f       id of the face to consider
 * \param[in]       center  point used for centering the set of basis functions
 * \param[in, out]  cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_mono_face_basis_setup(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f,
                       const cs_real_t          center[3],
                       cs_cell_builder_t       *cb)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ |
                       CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DIAM));

  for (int k = 0; k < 3; k++)
    bf->center[k] = center[k];

  /* Setting two generators of the face. They will be used as local axes.
   * Consider computing inertial axes and choosing them.
   * We chose the most "orthogonal" axes (this proved to lead to satisfactory
   * results), that is the two that give the minimum dot product.
   */

  cs_real_3_t *xexf_vect = cb->vectors;
  cs_real_t  tmp_meas = 0.;
  short int idx = 0;
  for (int j = cm->f2e_idx[f]; j < cm->f2e_idx[f+1]; j++)
    cs_math_3_length_unitv(center, cm->edge[cm->f2e_ids[j]].center,
                           &tmp_meas, xexf_vect[idx++]);

  cs_real_t min = 1.;
  short int e1 = 0, e2 = 1; // First selected couple
  const short int n_e = cm->f2e_idx[f+1] - cm->f2e_idx[f];
  for (short int i = 0; i < n_e; i++) {
    for (short int j = i+1; j < n_e; j++) {
      const cs_real_t  dp_ij = fabs(_dp3(xexf_vect[i], xexf_vect[j]));
      if (dp_ij < min)
        e1 = i, e2 = j, min = dp_ij;
    }
  }

  /* Assign axis */
  const cs_real_t  odiam = 1. / cm->f_diam[f];
  bf->axis[0].meas = bf->axis[1].meas = odiam;
  for (int k = 0; k < 3; k++) {
    bf->axis[0].unitv[k] = xexf_vect[e1][k];
    bf->axis[1].unitv[k] = xexf_vect[e2][k];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_BASIS_FUNC_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s : cell_id: %d, face: %d\n",
                __func__, cm->c_id, f);
  cs_log_printf(CS_LOG_DEFAULT,
                " Face axis             unitv                      meas\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n",
                bf->axis[0].unitv[0], bf->axis[0].unitv[1],
                bf->axis[0].unitv[2], bf->axis[0].meas,
                bf->axis[1].unitv[0], bf->axis[1].unitv[1],
                bf->axis[1].unitv[2], bf->axis[1].meas);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the face-basis function accordingly to the requested type
 *         Compute the inertia matrix related to the current face
 *
 * \param[in, out]  pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]       cm      pointer to a cs_cell_mesh_t
 * \param[in]       f       id of the face to consider
 * \param[in]       center  point used for centering the set of basis functions
 * \param[in, out]  cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_iner_face_basis_setup(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f,
                       const cs_real_t          center[3],
                       cs_cell_builder_t       *cb)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ));

  /* Initialization using the monomial basis */
  _mono_face_basis_setup(bf, cm, f, center, cb);

  /* Compute the covariance matrix related to the current face */
  cs_real_t  cov[3];
  cs_compute_face_covariance_tensor(cm, f, bf->axis[0], bf->axis[1], center,
                                    cov);

  /* Advanced parameters for controlling the algorithm */
  const double  tolerance = 1e-6;
  const cs_real_t cov_score = fabs(cov[1])/sqrt(cov[0]*cov[2]);

  if (cov_score < tolerance)
    /* axis is already an approximation of the inertial axes, that's ok
     * for us */
    return;

  /* Characteristic polynomial: l^2 - tr(A) l + det(A) */
  const cs_real_t  tr  = cov[0] + cov[2];
  const cs_real_t  det = cov[0]*cov[2] - cov[1]*cov[1];
  const cs_real_t  discrim = sqrt( tr*tr - 4.*det );

  /* First eigenvalue */
  const cs_real_t lambda = ( tr + discrim)*0.5;
  const cs_real_t mu0 = cov[0] - lambda, mu1 = cov[1];
  const double nrm2d = 1. / sqrt(mu0*mu0 + mu1*mu1);
  const double c0 = mu0 * nrm2d, c1 = mu1 * nrm2d;

  /* Assign axis */
  cs_real_3_t  c1a0_c0a1;
  for(int k = 0; k < 3; k++)
    c1a0_c0a1[k] = -c1 * bf->axis[0].unitv[k] +
                    c0 * bf->axis[1].unitv[k];

  cs_real_3_t  c0a0_c1a1;
  for(int k = 0; k < 3; k++)
    c0a0_c1a1[k] =  c0 * bf->axis[0].unitv[k] +
                    c1 * bf->axis[1].unitv[k];

  cs_nvec3(c1a0_c0a1, bf->axis);  // First axis
  cs_nvec3(c0a0_c1a1, bf->axis + 1);  // Second axis

  const cs_real_t  odiam = 1. / cm->f_diam[f];
  bf->axis[0].meas = odiam;
  bf->axis[1].meas = odiam;

#if defined(DEBUG) && !defined(NDEBUG) && CS_BASIS_FUNC_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s : cell_id: %d, face: %d\n",
                __func__, cm->c_id, f);
  cs_log_printf(CS_LOG_DEFAULT,
                " Face axis             unitv                      meas\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n"
                "           |% 5.3e % 5.3e % 5.3e|  % 5.3e\n",
                bf->axis[0].unitv[0], bf->axis[0].unitv[1],
                bf->axis[0].unitv[2], bf->axis[0].meas,
                bf->axis[1].unitv[0], bf->axis[1].unitv[1],
                bf->axis[1].unitv[2], bf->axis[1].meas);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to cell DoFs at a given point
 *         Case of a polynomial basis up to order 0
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_k0_eval_all_at_point(const void           *pbf,
                      const cs_real_t       coords[3],
                      cs_real_t            *eval)
{
  CS_UNUSED(coords);
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(eval != NULL);
  eval[0] = bf->phi0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to cell DoFs at a given point
 *         Case of a polynomial basis up to order 0
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_k0_eval_at_point(const void           *pbf,
                  const cs_real_t       coords[3],
                  short int             start,
                  short int             end,
                  cs_real_t            *eval)
{
  CS_UNUSED(coords);
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(eval != NULL);
  assert(start == 0 && end == 1);

  if (start >= end)
    return;

  eval[0] = bf->phi0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned cell basis functions in
 *         the case k=0 (constant value)
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      id      id of the element to consider
 */
/*----------------------------------------------------------------------------*/

static void
_ck0_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id)
{
  CS_UNUSED(id);
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(bf->size);

  cs_real_t  *pval = bf->projector->val;

  /* Only one entry in this case */
  bf->projector->n_rows = bf->projector->n_cols = 1;
  pval[0] = cm->vol_c * bf->phi0 * bf->phi0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned face basis functions in
 *         the case k=0 (constant value)
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]       f      id of the face in the cell numbering
 */
/*----------------------------------------------------------------------------*/

static void
_fk0_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(bf->size);

  cs_real_t  *pval = bf->projector->val;

  /* Only one entry in this case */
  bf->projector->n_rows = bf->projector->n_cols = 1;
  pval[0] = cm->face[f].meas * bf->phi0 * bf->phi0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the Modified Choloesky factorization of the projection
 *         matrix (mass matrix) related to the basis function
 *         Case of 0th order on cells or faces
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_k0_compute_facto(void       *pbf)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(bf->projector != NULL);
  assert(bf->projector->val[0] > 0);

  /* Mass matrix is a 1x1 matrix ! */
  if (bf->facto_max_size < 1) {
    BFT_REALLOC(bf->facto, 1, cs_real_t);
    bf->facto_max_size = 1;
  }
  bf->facto[0] = 1/bf->projector->val[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project an array on the polynomial basis functions of order 0.
 *         This results from the application of a Modified Choloesky
 *         factorization which should be performed before calling this function.
 *         The input array is defined as follows:
 *         int_elem v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

static void
_k0_project(const void              *pbf,
            const cs_real_t         *array,
            cs_real_t               *dof)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  assert(array != NULL && dof != NULL);

  dof[0] = bf->facto[0] * array[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to cell DoFs at a given point
 *         Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_ck1_eval_all_at_point(const void           *pbf,
                       const cs_real_t       coords[3],
                       cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  eval[0] = bf->phi0;
  eval[1] = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  eval[2] = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;
  eval[3] = _dp3(r, bf->axis[2].unitv) * bf->axis[2].meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a set of basis functions related to cell DoFs at a given
 *         point. Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_ck1_eval_at_point(const void           *pbf,
                   const cs_real_t       coords[3],
                   short int             start,
                   short int             end,
                   cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  assert(start > -1 && end < 5);
  short int shift = 0;
  for (short int i = start; i < end; i++, shift++) {
    if (i == 0)
      eval[shift] = bf->phi0;
    else
      eval[shift] = _dp3(r, bf->axis[i-1].unitv) * bf->axis[i-1].meas;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned cell basis functions in
 *         the case k=1 (up to 1st order polynomial basis)
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      id      id of the element to consider
 */
/*----------------------------------------------------------------------------*/

static void
_ck1_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id)
{
  CS_UNUSED(id);
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  const int n_rows = CK1_SIZE;
  const int n_gpts = 4;

  cs_real_t  phi_eval[CK1_SIZE], weights[4];
  cs_real_3_t  gpts[4];
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(bf->size);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

  /* First row (or first column since projector is symmetric) is easy.
     So, 1st row is computed differently */
#if defined(DEBUG) && !defined(NDEBUG)
  bf->eval_all_at_point(bf, cm->xc, phi_eval);
  pval[0] = cm->vol_c * phi_eval[0];
  pval[1] = cm->vol_c * phi_eval[1];
  pval[2] = cm->vol_c * phi_eval[2];
  pval[3] = cm->vol_c * phi_eval[3];
#else
  pval[0] = cm->vol_c;
#endif

  switch (cm->type) {

  case FVM_CELL_TETRA:
    cs_quadrature_tet_4pts(cm->xv, cm->xv +3, cm->xv +6, cm->xv +9, cm->vol_c,
                           gpts, weights);
    _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                 1, n_rows, phi_eval, pval);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:

    /* Loop on faces */
    for (short int f = 0; f < cm->n_fc; f++) {

      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_ef = end - start; /* #vertices (= #edges) */
      const short int *f2e_ids = cm->f2e_ids + start;
      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_1ov3 * cm->hfc[f];

      switch (n_ef) { /* Optimized version for triangles */
      case CS_TRIANGLE_CASE: /* Triangle */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          cs_quadrature_tet_4pts(cm->xv + 3*v0, cm->xv + 3*v1,
                                 cm->xv + 3*v2, cm->xc, hf_coef * pfq.meas,
                                 gpts, weights);
          _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                       1, n_rows, phi_eval, pval);
        }
        break;

      default:
        {
          assert(n_ef > 3);
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const short int v0 = cm->e2v_ids[2*e0];
            const short int v1 = cm->e2v_ids[2*e0+1];

            cs_quadrature_tet_4pts(cm->xv + 3*v0, cm->xv + 3*v1,
                                   pfq.center, cm->xc, hf_coef * tef[e],
                                   gpts, weights);
            _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                         1, n_rows, phi_eval, pval);

          }
        }
        break;

      } /* End of switch on n_ef */
    } /* Loop on cell faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Unknown cell-type.\n", __func__);
    break;
  }

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the Modified Choloesky factorization of the projection
 *         matrix (mass matrix) related to the basis function
 *         Case of 2nd order on cells
 *
 * \param[in, out]  pbf    pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ck1_compute_facto(void       *pbf)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(bf->projector != NULL);

  /* Mass matrix is a 4x4 matrix */
  if (bf->facto_max_size < 10) {
    bf->facto_max_size = 10;
    BFT_REALLOC(bf->facto, 10, cs_real_t);
  }
  cs_sdm_44_ldlt_compute(bf->projector, bf->facto);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project an array on the polynomial cell basis functions up to first
 *         order. This results from the application of a Modified Choloesky
 *         factorization which should be performed before calling this function.
 *         The input array is defined as follows:
 *         int_elem v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

static inline void
_ck1_project(const void              *pbf,
             const cs_real_t         *array,
             cs_real_t               *dof)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  assert(array != NULL && dof != NULL);

  /* array stands for the rhs and dof the solution */
  cs_sdm_44_ldlt_solve(bf->facto, array, dof);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned cell basis functions in
 *         the case k=2 (up to 1st order polynomial basis)
 *         1, x, y, z, xx, xy, xz, yy, yz, zz
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      id      id of the element to consider
 */
/*----------------------------------------------------------------------------*/

static void
_ck2_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id)
{
  CS_UNUSED(id);
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  const int n_rows = CK2_SIZE;
  const int n_gpts = 15;

  cs_real_t  phi_eval[CK2_SIZE], weights[15];
  cs_real_3_t  gpts[15];
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(bf->size);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

  switch (cm->type) {

  case FVM_CELL_TETRA:
    cs_quadrature_tet_15pts(cm->xv, cm->xv +3, cm->xv +6, cm->xv +9, cm->vol_c,
                            gpts, weights);
    _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                 0, n_rows, phi_eval, pval);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:

    /* Loop on faces */
    for (short int f = 0; f < cm->n_fc; f++) {

      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_ef = end - start; /* #vertices (= #edges) */
      const short int *f2e_ids = cm->f2e_ids + start;
      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_1ov3 * cm->hfc[f];

      switch (n_ef) { /* Optimized version for triangles */
      case CS_TRIANGLE_CASE: /* Triangle */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          cs_quadrature_tet_15pts(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2,
                                  cm->xc, hf_coef * pfq.meas,
                                  gpts, weights);
          _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                       0, n_rows, phi_eval, pval);
        }
        break;

      default:
        {
          assert(n_ef > 3);
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const short int v0 = cm->e2v_ids[2*e0];
            const short int v1 = cm->e2v_ids[2*e0+1];

            cs_quadrature_tet_15pts(cm->xv + 3*v0, cm->xv + 3*v1,
                                    pfq.center, cm->xc, hf_coef * tef[e],
                                    gpts, weights);
            _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                         0, n_rows, phi_eval, pval);

          }
        }
        break;

      } /* End of switch on n_ef */
    } /* Loop on cell faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Unknown cell-type.\n", __func__);
    break;
  }

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to cell DoFs at a given point
 *         Case of a polynomial basis for arbitrary order
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cka_eval_all_at_point(const void           *pbf,
                       const cs_real_t       coords[3],
                       cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  eval[0] = bf->phi0;
  eval[1] = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  eval[2] = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;
  eval[3] = _dp3(r, bf->axis[2].unitv) * bf->axis[2].meas;

  cs_real_t  *eval_ho = eval + 4;
  for (int i = 0; i < bf->n_deg_elts; i++)
    eval_ho[i] = _int_pow(eval[1], bf->deg[3*i  ]) *
                 _int_pow(eval[2], bf->deg[3*i+1]) *
                 _int_pow(eval[3], bf->deg[3*i+2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a set of basis functions related to cell DoFs at a given
 *         point. Case of a polynomial basis for arbitrary order
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cka_eval_at_point(const void           *pbf,
                   const cs_real_t       coords[3],
                   short int             start,
                   short int             end,
                   cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);
  assert(start > -1 && end < 4 + bf->n_deg_elts + 1);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  cs_real_t  p1_eval[3] = {0, 0, 0};
  for (short int i = 0; i < 3; i++)
    p1_eval[i] = _dp3(r, bf->axis[i].unitv) * bf->axis[i].meas;

  short int shift = 0;
  for (short int i = start; i < end; i++, shift++) {
    if (i == 0)
      eval[shift] = bf->phi0;
    else if (i < 4)
      eval[shift] = p1_eval[i-1];
    else {
      short int j = i - 4;
      eval[shift] = _int_pow(p1_eval[0], bf->deg[3*j  ]) *
                    _int_pow(p1_eval[1], bf->deg[3*j+1]) *
                    _int_pow(p1_eval[2], bf->deg[3*j+2]);
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned cell basis functions in
 *         the case k=2 (up to 1st order polynomial basis)
 *         1, x, y, z, xx, xy, xz, yy, yz, zz
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      id      id of the element to consider
 */
/*----------------------------------------------------------------------------*/

static void
_cka_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          id)
{
  CS_UNUSED(id);
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  const int n_rows = bf->size;
  /* Max quadrature available 15pts. Problem if k > 2 */
  const int n_gpts = bf->n_gpts_tetra;

  cs_real_t  *phi_eval = NULL, *weights = NULL;
  cs_real_3_t  *gpts = NULL;
  cs_real_t  _phi_eval[CKA_SIZE], _weights[CKA_SIZE];
  cs_real_3_t  _gpts[CKA_SIZE];

  if (n_rows > CKA_SIZE)
    BFT_MALLOC(phi_eval, n_rows, cs_real_t);
  else
    phi_eval = _phi_eval;
  if (n_gpts > CKA_SIZE) {
    BFT_MALLOC(weights, n_gpts, cs_real_t);
    BFT_MALLOC(gpts, n_gpts, cs_real_3_t);
  }
  else {
    weights = _weights;
    gpts = _gpts;
  }

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(n_rows);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

  switch (cm->type) {

  case FVM_CELL_TETRA:
    bf->quadrature_tetra(cm->xv, cm->xv +3, cm->xv +6, cm->xv +9, cm->vol_c,
                         gpts, weights);
    _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                 0, n_rows, phi_eval, pval);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:

    /* Loop on faces */
    for (short int f = 0; f < cm->n_fc; f++) {

      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_ef = end - start; /* #vertices (= #edges) */
      const short int *f2e_ids = cm->f2e_ids + start;
      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_1ov3 * cm->hfc[f];

      switch (n_ef) { /* Optimized version for triangles */
      case CS_TRIANGLE_CASE: /* Triangle */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          bf->quadrature_tetra(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2,
                               cm->xc, hf_coef * pfq.meas,
                               gpts, weights);
          _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                       0, n_rows, phi_eval, pval);
        }
        break;

      default:
        {
          assert(n_ef > 3);
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const short int v0 = cm->e2v_ids[2*e0];
            const short int v1 = cm->e2v_ids[2*e0+1];

            bf->quadrature_tetra(cm->xv + 3*v0, cm->xv + 3*v1,
                                 pfq.center, cm->xc, hf_coef * tef[e],
                                 gpts, weights);
            _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                         0, n_rows, phi_eval, pval);

          }
        }
        break;

      } /* End of switch on n_ef */
    } /* Loop on cell faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Unknown cell-type.\n", __func__);
    break;
  }

  /* Free allocated buffers */
  if (n_rows > CKA_SIZE)
    BFT_FREE(phi_eval);
  if (n_gpts > CKA_SIZE) {
    BFT_FREE(weights);
    BFT_FREE(gpts);
  }

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the Modified Choloesky factorization of the projection
 *         matrix (mass matrix) related to the basis function
 *         Case of arbitrary order on cells/faces
 *
 * \param[in, out]  pbf    pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_ka_compute_facto(void       *pbf)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(bf->projector != NULL);

  /* Mass matrix is a nxn matrix. Facto is allocated to a size which is
     greater than the one requested in order to store a temporary buffer
     used for building the factorization */
  const int facto_size = ((bf->size + 1)*bf->size)/2;
  if (bf->facto_max_size < facto_size + bf->size) {
    bf->facto_max_size = facto_size + bf->size;
    BFT_REALLOC(bf->facto, facto_size + bf->size, cs_real_t);
  }
  cs_sdm_ldlt_compute(bf->projector, bf->facto, bf->facto + facto_size);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project an array on the polynomial cell/face basis functions up to
 *         arbitrary order. This results from the application of a Modified
 *         Choloesky factorization which should be performed before calling
 *         this function.
 *         The input array is defined as follows:
 *         int_elem v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

static inline void
_ka_project(const void              *pbf,
            const cs_real_t         *array,
            cs_real_t               *dof)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  assert(array != NULL && dof != NULL);

  /* array stands for the rhs and dof the solution */
  cs_sdm_ldlt_solve(bf->size, bf->facto, array, dof);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the projection matrix related to an arbitrary order polynomial
 *         cell basis functions
 *
 * \param[in]    pbf     pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_cka_dump_projector(const void      *pbf)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  const cs_real_t  *pval = bf->projector->val;
  const char  *labels[20] = {"0", "x", "y", "z",
                             "x.x", "x.y", "x.z", "y.y", "y.z", "z.z",
                             "x3" , "x2y", "x2z", "xy2", "xyz", "xz2", "y3",
                             "y2z", "yz2", "z3"};

  if (bf->size == 0)
    return;

  cs_log_printf(CS_LOG_DEFAULT, "%6s %11s", " ", labels[0]);
  if (bf->size > 20) {

    for (int i = 1; i < 20; i++)
      cs_log_printf(CS_LOG_DEFAULT, " %11s", labels[i]);
    for (int i = 0; i < bf->size; i++) {
      if (i < 20)
        cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", labels[i]);
      else
        cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", " ");
      for (int j = 0; j < bf->size; j++)
        cs_log_printf(CS_LOG_DEFAULT, " % .4e", pval[j]);
      pval += bf->size;
    }

  }
  else {

    for (int i = 1; i < bf->size; i++)
      cs_log_printf(CS_LOG_DEFAULT, " %11s", labels[i]);
    for (int i = 0; i < bf->size; i++) {
      cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", labels[i]);
      for (int j = 0; j < bf->size; j++)
        cs_log_printf(CS_LOG_DEFAULT, " % .4e", pval[j]);
      pval += bf->size;
    }

  }
  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all gradient of the basis functions related to cell DoFs
 *         at a given point
 *         Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cgk1_eval_all_at_point(const void           *pbf,
                        const cs_real_t       coords[3],
                        cs_real_t            *eval)
{
  CS_UNUSED(coords);
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(eval != NULL);
  /* Gradient of constant is zero */
  eval[0] = eval[1] = eval[2] = 0;

  /* Gradient of affine functions is a constant vector */
  eval[3] = bf->axis[0].unitv[0] * bf->axis[0].meas;
  eval[4] = bf->axis[0].unitv[1] * bf->axis[0].meas;
  eval[5] = bf->axis[0].unitv[2] * bf->axis[0].meas;

  eval[6] = bf->axis[1].unitv[0] * bf->axis[1].meas;
  eval[7] = bf->axis[1].unitv[1] * bf->axis[1].meas;
  eval[8] = bf->axis[1].unitv[2] * bf->axis[1].meas;

  eval[ 9] = bf->axis[2].unitv[0] * bf->axis[2].meas;
  eval[10] = bf->axis[2].unitv[1] * bf->axis[2].meas;
  eval[11] = bf->axis[2].unitv[2] * bf->axis[2].meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the gradient for a set of basis functions related to cell
 *         DoFs at a given point.
 *         Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cgk1_eval_at_point(const void           *pbf,
                    const cs_real_t       coords[3],
                    short int             start,
                    short int             end,
                    cs_real_t            *eval)
{
  CS_UNUSED(coords);
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(eval != NULL);
  assert(start > -1 && end < 5);

  short int shift = 0;
  for (short int i = start; i < end; i++, shift += 3) {

    if (i == 0) {
      /* Gradient of constant is zero */
      eval[shift] = eval[shift+1] = eval[shift+2] = 0;
    }
    else {
      /* Gradient of affine functions is a constant vector */
      eval[shift  ] = bf->axis[i-1].unitv[0] * bf->axis[i-1].meas;
      eval[shift+1] = bf->axis[i-1].unitv[1] * bf->axis[i-1].meas;
      eval[shift+2] = bf->axis[i-1].unitv[2] * bf->axis[i-1].meas;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all gradient of the basis functions related to cell DoFs
 *         at a given point
 *         Case of a polynomial basis with order greater than 2
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cgka_eval_all_at_point(const void           *pbf,
                        const cs_real_t       coords[3],
                        cs_real_t            *eval)
{
  assert(coords != NULL && eval != NULL);

  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  /* Gradient of constant is zero */
  eval[0] = eval[1] = eval[2] = 0;

  /* Gradient of affine functions is a constant vector corresponding to
     to the inertial axis */
  eval[3] = bf->axis[0].unitv[0] * bf->axis[0].meas;
  eval[4] = bf->axis[0].unitv[1] * bf->axis[0].meas;
  eval[5] = bf->axis[0].unitv[2] * bf->axis[0].meas;

  eval[6] = bf->axis[1].unitv[0] * bf->axis[1].meas;
  eval[7] = bf->axis[1].unitv[1] * bf->axis[1].meas;
  eval[8] = bf->axis[1].unitv[2] * bf->axis[1].meas;

  eval[ 9] = bf->axis[2].unitv[0] * bf->axis[2].meas;
  eval[10] = bf->axis[2].unitv[1] * bf->axis[2].meas;
  eval[11] = bf->axis[2].unitv[2] * bf->axis[2].meas;

  /* Set higher order terms */
  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  /* Monomial for x, y and z */
  const cs_real_t  mx = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  const cs_real_t  my = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;
  const cs_real_t  mz = _dp3(r, bf->axis[2].unitv) * bf->axis[2].meas;

  cs_real_t  *eval_ho = eval + 12;
  for (int i = 0; i < bf->n_deg_elts; i++) {

    const short int  x_expo = bf->deg[3*i];
    const short int  y_expo = bf->deg[3*i+1];
    const short int  z_expo = bf->deg[3*i+2];

    cs_real_t  *grd = eval_ho + 3*i;
    grd[0] = grd[1] = grd[2] = 0;

    if (x_expo > 0) {
      const cs_real_t  phi = x_expo * bf->axis[0].meas *
        _int_pow(mx, x_expo-1) * _int_pow(my, y_expo) * _int_pow(mz, z_expo);
      grd[0] += phi * bf->axis[0].unitv[0];
      grd[1] += phi * bf->axis[0].unitv[1];
      grd[2] += phi * bf->axis[0].unitv[2];
    }
    if (y_expo > 0) {
      const cs_real_t  phi = y_expo * bf->axis[1].meas *
        _int_pow(mx, x_expo) * _int_pow(my, y_expo-1) * _int_pow(mz, z_expo);
      grd[0] += phi * bf->axis[1].unitv[0];
      grd[1] += phi * bf->axis[1].unitv[1];
      grd[2] += phi * bf->axis[1].unitv[2];
    }
    if (z_expo > 0) {
      const cs_real_t  phi = z_expo * bf->axis[2].meas *
        _int_pow(mx, x_expo) * _int_pow(my, y_expo) * _int_pow(mz, z_expo-1);
      grd[0] += phi * bf->axis[2].unitv[0];
      grd[1] += phi * bf->axis[2].unitv[1];
      grd[2] += phi * bf->axis[2].unitv[2];
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the gradient of a set of basis functions related to cell
 *         DoFs at a given point.
 *         Case of a polynomial basis with order greater than 2.
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_cgka_eval_at_point(const void           *pbf,
                    const cs_real_t       coords[3],
                    short int             start,
                    short int             end,
                    cs_real_t            *eval)
{
  assert(coords != NULL && eval != NULL);
  assert(start > -1);

  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  /* Set higher order terms */
  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  /* Monomial for x, y and z */
  const cs_real_t  mx = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  const cs_real_t  my = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;
  const cs_real_t  mz = _dp3(r, bf->axis[2].unitv) * bf->axis[2].meas;

  short int shift = 0;
  for (short int i = start; i < end; i++, shift += 3) {

    if (i == 0) {
      /* Gradient of constant is zero */
      eval[shift] = eval[shift+1] = eval[shift+2] = 0;
    }
    else if (i < 4) {
      /* Gradient of affine functions is a constant vector */
      eval[shift  ] = bf->axis[i-1].unitv[0] * bf->axis[i-1].meas;
      eval[shift+1] = bf->axis[i-1].unitv[1] * bf->axis[i-1].meas;
      eval[shift+2] = bf->axis[i-1].unitv[2] * bf->axis[i-1].meas;
    }
    else {

      const int j = i -4; /* shift in n_deg_elts */
      assert(j > -1);
      const short int  x_expo = bf->deg[3*j];
      const short int  y_expo = bf->deg[3*j+1];
      const short int  z_expo = bf->deg[3*j+2];

      /* Set to zero before proceeding to an accumulation */
      eval[shift] = eval[shift+1] = eval[shift+2] = 0;

      if (x_expo > 0) {
        const cs_real_t  phi = x_expo * bf->axis[0].meas *
          _int_pow(mx, x_expo-1) * _int_pow(my, y_expo) * _int_pow(mz, z_expo);
        eval[shift  ] += phi * bf->axis[0].unitv[0];
        eval[shift+1] += phi * bf->axis[0].unitv[1];
        eval[shift+2] += phi * bf->axis[0].unitv[2];
      }
      if (y_expo > 0) {
        const cs_real_t  phi = y_expo * bf->axis[1].meas *
          _int_pow(mx, x_expo) * _int_pow(my, y_expo-1) * _int_pow(mz, z_expo);
        eval[shift  ] += phi * bf->axis[1].unitv[0];
        eval[shift+1] += phi * bf->axis[1].unitv[1];
        eval[shift+2] += phi * bf->axis[1].unitv[2];
      }
      if (z_expo > 0) {
        const cs_real_t  phi = z_expo * bf->axis[2].meas *
          _int_pow(mx, x_expo) * _int_pow(my, y_expo) * _int_pow(mz, z_expo-1);
        eval[shift  ] += phi * bf->axis[2].unitv[0];
        eval[shift+1] += phi * bf->axis[2].unitv[1];
        eval[shift+2] += phi * bf->axis[2].unitv[2];
      }

    } /* Higher order terms */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to face DoFs at a given point
 *         Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_fk1_eval_all_at_point(const void           *pbf,
                       const cs_real_t       coords[3],
                       cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  eval[0] = bf->phi0;
  eval[1] = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  eval[2] = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to face DoFs at a given point
 *         Case of a polynomial basis up to order 1
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_fk1_eval_at_point(const void           *pbf,
                   const cs_real_t       coords[3],
                   short int             start,
                   short int             end,
                   cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);
  assert(start > -1 && end < 4);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  short int shift = 0;
  for (short int i = start; i < end; i++, shift++) {
    if (i == 0)
      eval[shift] = bf->phi0;
    else
      eval[shift] = _dp3(r, bf->axis[i-1].unitv) * bf->axis[i-1].meas;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned face basis functions in
 *         the case k=1 (up to 1st order polynomial basis)
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      f       id of the face to consider (in cell numbering)
 */
/*----------------------------------------------------------------------------*/

static void
_fk1_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  /* First row (or first column since projector is symmetric) is easy.
     So, 1st row is computed differently */
  const int n_rows = FK1_SIZE;
  const int n_gpts = 3;
  const cs_quant_t  pfq = cm->face[f];

  cs_real_t  phi_eval[FK1_SIZE], weights[3];
  cs_real_3_t  gpts[3];

  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(n_rows);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

#if defined(DEBUG) && !defined(NDEBUG)
  bf->eval_all_at_point(bf, pfq.center, phi_eval);
  pval[0] = pfq.meas * phi_eval[0];
  pval[1] = pfq.meas * phi_eval[1];
  pval[2] = pfq.meas * phi_eval[2];
#else
  pval[0] = pfq.meas * bf->phi0;;
#endif

  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int n_ef = end - start; /* #vertices (= #edges) */
  const short int *f2e_ids = cm->f2e_ids + start;

  switch (n_ef) { /* Optimized version for triangles */

  case CS_TRIANGLE_CASE: /* Triangle */
    {
      short int  v0, v1, v2;
      cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

      cs_quadrature_tria_3pts(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2,
                              pfq.meas, gpts, weights);
      weights[2] = weights[1] = weights[0];
      _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                   1, n_rows, phi_eval, pval);
    }
    break;

  default:
    {
      assert(n_ef > 3);
      const double  *tef = cm->tef + start;

      for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

        // Edge-related variables
        const short int e0  = f2e_ids[e];
        const short int v0 = cm->e2v_ids[2*e0];
        const short int v1 = cm->e2v_ids[2*e0+1];

        cs_quadrature_tria_3pts(cm->xv + 3*v0, cm->xv + 3*v1, pfq.center,
                                tef[e], gpts, weights);
        weights[2] = weights[1] = weights[0];
        _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                     1, n_rows, phi_eval, pval);

      }
    }
    break;

  } /* End of switch on n_ef */

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the Modified Choloesky factorization of the projection
 *         matrix (mass matrix) related to the basis function
 *         Case of 1st order on faces
 *
 * \param[in, out]  pbf    pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_fk1_compute_facto(void       *pbf)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(bf->projector != NULL);

  /* Mass matrix is a 3x3 matrix. Factorization is stored with 6 values */
  if (bf->facto_max_size < 6) {
    bf->facto_max_size = 6;
    BFT_REALLOC(bf->facto, 6, cs_real_t);
  }
  cs_sdm_33_ldlt_compute(bf->projector, bf->facto);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project an array on the polynomial face basis functions up to first
 *         order. This results from the application of a Modified Choloesky
 *         factorization which should be performed before calling this function.
 *         The input array is defined as follows:
 *         int_face v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

static inline void
_fk1_project(const void              *pbf,
             const cs_real_t         *array,
             cs_real_t               *dof)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  assert(array != NULL && dof != NULL);

  /* array stands for the rhs and dof the solution */
  cs_sdm_33_ldlt_solve(bf->facto, array, dof);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned face basis functions in
 *         the case k=2 (up to 2nd order polynomial basis)
 *         1, xf, yf, xf.xf, xf.yf, yf.yf
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      f       id of the face to consider (in cell numbering)
 */
/*----------------------------------------------------------------------------*/

static void
_fk2_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  /* First row (or first column since projector is symmetric) is easy.
     So, 1st row is computed differently */
  const int n_rows = FK2_SIZE;
  const int n_gpts = 7;
  const cs_quant_t  pfq = cm->face[f];

  cs_real_t  phi_eval[FK2_SIZE], weights[7];
  cs_real_3_t  gpts[7];

  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(n_rows);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int n_ef = end - start; /* #vertices (= #edges) */
  const short int *f2e_ids = cm->f2e_ids + start;

  switch (n_ef) { /* Optimized version for triangles */

  case CS_TRIANGLE_CASE: /* Triangle */
    {
      short int  v0, v1, v2;
      cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

      cs_quadrature_tria_7pts(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2,
                              pfq.meas, gpts, weights);
      _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                   0, n_rows, phi_eval, pval);
    }
    break;

  default:
    {
      assert(n_ef > 3);
      const double  *tef = cm->tef + start;

      for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

        // Edge-related variables
        const short int e0  = f2e_ids[e];
        const short int v0 = cm->e2v_ids[2*e0];
        const short int v1 = cm->e2v_ids[2*e0+1];

        cs_quadrature_tria_7pts(cm->xv + 3*v0, cm->xv + 3*v1, pfq.center,
                                tef[e], gpts, weights);
        _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                     0, n_rows, phi_eval, pval);

      }
    }
    break;

  } /* End of switch on n_ef */

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the Modified Choloesky factorization of the projection
 *         matrix (mass matrix) related to the basis function
 *         Case of 2nd order on faces
 *
 * \param[in, out]  pbf    pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_fk2_compute_facto(void       *pbf)
{
  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;
  assert(bf->projector != NULL);

  /* Mass matrix is a 6x6 matrix */
  if (bf->facto_max_size < 21) {
    bf->facto_max_size = 21;
    BFT_REALLOC(bf->facto, 21, cs_real_t);
  }

  cs_sdm_66_ldlt_compute(bf->projector, bf->facto);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project an array on the polynomial face basis functions up to second
 *         order. This results from the application of a Modified Choloesky
 *         factorization which should be performed before calling this function.
 *         The input array is defined as follows:
 *         int_face v.phi_i for all i in the basis. v is a function to estimate
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      array   array to project
 * \param[in, out] dof     projection of the array (= DoF values in this basis)
 */
/*----------------------------------------------------------------------------*/

static inline void
_fk2_project(const void              *pbf,
             const cs_real_t         *array,
             cs_real_t               *dof)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  assert(array != NULL && dof != NULL);

  /* array stands for the rhs and dof the solution */
  cs_sdm_66_ldlt_solve(bf->facto, array, dof);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate all basis functions related to face DoFs at a given point
 *         Case of a polynomial basis up to arbitrary order.
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_fka_eval_all_at_point(const void           *pbf,
                       const cs_real_t       coords[3],
                       cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  eval[0] = bf->phi0;
  eval[1] = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  eval[2] = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;

  cs_real_t  *eval_ho = eval + 3;
  for (int i = 0; i < bf->n_deg_elts; i++)
    eval_ho[i] =
      _int_pow(eval[1], bf->deg[2*i]) * _int_pow(eval[2], bf->deg[2*i+1]);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a set of basis functions related to face DoFs at a given
 *         point. Case of a polynomial basis up to arbitrary order.
 *
 * \param[in]      pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      coords  point coordinates
 * \param[in]      start   starts evaluating basis function from this id
 * \param[in]      end     ends evaluating basis function at this id-1
 * \param[in, out] eval    vector containing the evaluations of the functions
 */
/*----------------------------------------------------------------------------*/

static inline void
_fka_eval_at_point(const void           *pbf,
                   const cs_real_t       coords[3],
                   short int             start,
                   short int             end,
                   cs_real_t            *eval)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;

  assert(coords != NULL && eval != NULL);
  assert(start > -1 && end < bf->n_deg_elts + 1);

  const cs_real_3_t  r = {coords[0] - bf->center[0],
                          coords[1] - bf->center[1],
                          coords[2] - bf->center[2]};

  cs_real_t  mono[2];
  mono[0] = _dp3(r, bf->axis[0].unitv) * bf->axis[0].meas;
  mono[1] = _dp3(r, bf->axis[1].unitv) * bf->axis[1].meas;

  short int  shift = 0;
  for (short int i = start; i < end; i++, shift++) {
    if (i == 0)
      eval[shift] = bf->phi0;
    else if (i < 3)
      eval[shift] = mono[i-1];
    else {
      const short int j = i - 3;
      eval[shift] =
        _int_pow(mono[0], bf->deg[2*j]) * _int_pow(mono[1], bf->deg[2*j+1]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projector to the space spanned face basis functions in
 *         the case k=2 (up to 2nd order polynomial basis)
 *         1, xf, yf, xf.xf, xf.yf, yf.yf
 *
 * \param[in, out] pbf     pointer to be cast into a cs_basis_func_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t
 * \param[in]      f       id of the face to consider (in cell numbering)
 */
/*----------------------------------------------------------------------------*/

static void
_fka_compute_projector(void                    *pbf,
                       const cs_cell_mesh_t    *cm,
                       const short int          f)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  cs_basis_func_t  *bf = (cs_basis_func_t *)pbf;

  const cs_quant_t  pfq = cm->face[f];
  const int n_rows = bf->size;
  /* Max quadrature available 7pts. Problem if k > 3 */
  const int n_gpts = bf->n_gpts_tria;

  cs_real_t  *phi_eval = NULL, *weights = NULL;
  cs_real_3_t  *gpts = NULL;
  cs_real_t  _phi_eval[FKA_SIZE], _weights[FKA_SIZE];
  cs_real_3_t  _gpts[FKA_SIZE];

  if (n_rows > FKA_SIZE)
    BFT_MALLOC(phi_eval, n_rows, cs_real_t);
  else
    phi_eval = _phi_eval;
  if (n_gpts > FKA_SIZE) {
    BFT_MALLOC(weights, n_gpts, cs_real_t);
    BFT_MALLOC(gpts, n_gpts, cs_real_3_t);
  }
  else {
    weights = _weights;
    gpts = _gpts;
  }

  if (bf->projector == NULL)
    bf->projector = cs_sdm_square_create(n_rows);

  /* Reset values */
  cs_sdm_square_init(n_rows, bf->projector);

  cs_real_t  *pval = bf->projector->val;

  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int n_ef = end - start; /* #vertices (= #edges) */
  const short int *f2e_ids = cm->f2e_ids + start;

  switch (n_ef) { /* Optimized version for triangles */

  case CS_TRIANGLE_CASE: /* Triangle */
    {
      short int  v0, v1, v2;
      cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

      bf->quadrature_tria(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2,
                          pfq.meas, gpts, weights);
      _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                   0, n_rows, phi_eval, pval);
    }
    break;

  default:
    {
      assert(n_ef > 3);
      const double  *tef = cm->tef + start;

      for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

        // Edge-related variables
        const short int e0  = f2e_ids[e];
        const short int v0 = cm->e2v_ids[2*e0];
        const short int v1 = cm->e2v_ids[2*e0+1];

        bf->quadrature_tria(cm->xv + 3*v0, cm->xv + 3*v1, pfq.center,
                            tef[e], gpts, weights);
        _add_contrib(n_gpts, (const cs_real_t (*)[3])gpts, weights, bf,
                     0, n_rows, phi_eval, pval);

      }
    }
    break;

  } /* End of switch on n_ef */

  /* Free allocated buffers */
  if (n_rows > FKA_SIZE)
    BFT_FREE(phi_eval);
  if (n_gpts > FKA_SIZE) {
    BFT_FREE(weights);
    BFT_FREE(gpts);
  }

  /* Projection matrix is symmetric by construction */
  _symmetrize_and_clean(n_rows, pval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the projection matrix related to 2nd order polynomial face
 *         basis functions
 *
 * \param[in]    pbf     pointer to be cast into a cs_basis_func_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_fka_dump_projector(const void      *pbf)
{
  const cs_basis_func_t  *bf = (const cs_basis_func_t *)pbf;
  const cs_real_t  *pval = bf->projector->val;
  const char  *labels[15] = {"0" , "x"  , "y"   ,
                             "x2", "xy" , "y2"  ,
                             "x3", "x2y", "xy2" , "y3",
                             "x4", "x3y", "x2y2", "xy3", "y4"};

  cs_log_printf(CS_LOG_DEFAULT, "%6s %11s", " ", labels[0]);
  if (bf->size > 15) {

    for (int i = 1; i < 15; i++)
      cs_log_printf(CS_LOG_DEFAULT, " %11s", labels[i]);
    for (int i = 0; i < bf->size; i++) {
      if (i < 15)
        cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", labels[i]);
      else
        cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", " ");
      for (int j = 0; j < bf->size; j++)
        cs_log_printf(CS_LOG_DEFAULT, " % .4e", pval[j]);
      pval += bf->size;
    }

  }
  else {

    for (int i = 1; i < bf->size; i++)
      cs_log_printf(CS_LOG_DEFAULT, " %11s", labels[i]);
    for (int i = 0; i < bf->size; i++) {
      cs_log_printf(CS_LOG_DEFAULT, "\n %6s ", labels[i]);
      for (int j = 0; j < bf->size; j++)
        cs_log_printf(CS_LOG_DEFAULT, " % .4e", pval[j]);
      pval += bf->size;
    }

  }
  cs_log_printf(CS_LOG_DEFAULT, "\n");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_basis_func_t structure
 *
 * \param[in]  flag     metadata related to the way of building basis functions
 * \param[in]  k        polynomial order
 * \param[in]  dim      2 or 3 w.r.t. the geometrical dimension
 *
 * \return a pointer to the new cs_basis_func_t
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_create(cs_flag_t      flag,
                     short int      k,
                     short int      dim)
{
  /* Structure to handle polynomial basis functions */
  cs_basis_func_t  *pbf = NULL;

  BFT_MALLOC(pbf, 1, cs_basis_func_t);

  pbf->flag = flag;
  pbf->poly_order = k;
  pbf->dim = dim;
  pbf->size = cs_math_binom(k + dim, dim);
  pbf->phi0 = 1;
  BFT_MALLOC(pbf->axis, dim, cs_nvec3_t);

  pbf->n_deg_elts = 0;
  pbf->deg = NULL;
  if (k > 1) {

    /* Remove the 3 first trivial basis functions in 2D
       and the 4 first trivial basis functions in 3D */
    pbf->n_deg_elts = pbf->size -(dim + 1);
    BFT_MALLOC(pbf->deg, pbf->n_deg_elts*dim, short int);

    if (dim == 3) {

      short int count = 0;
      for (short int ik = 2; ik < k+1; ik++) {
        for (short int e1 = ik; e1 > -1; e1--) {
          for (short int e2 = ik-e1; e2 > -1; e2--) {
            pbf->deg[count*dim    ] = e1;
            pbf->deg[count*dim + 1] = e2;
            pbf->deg[count*dim + 2] = ik - e1 -e2;
            count++;
          }
        }
      }
      assert(count == pbf->n_deg_elts);

    }
    else {

      assert(dim == 2);
      short int count = 0;
      for (short int ik = 2; ik < k+1; ik++) {
        for (short int e1 = ik; e1 > -1; e1--) {
          pbf->deg[count*dim    ] = e1;
          pbf->deg[count*dim + 1] = ik - e1;
          count++;
        }
      }
      assert(count == pbf->n_deg_elts);

    }

  } /* Scheme order > 1 */

  /* Handle projector (mass matrix) */
  pbf->projector = NULL;
  pbf->dump_projector = NULL;
  pbf->compute_projector = NULL;
  pbf->compute_factorization = NULL;
  pbf->facto_max_size = 0;
  pbf->facto = NULL;  /* array storing the matrix factorization */
  pbf->project = NULL;

  pbf->n_gpts_tria = 0;
  pbf->quadrature_tria = NULL;
  pbf->n_gpts_tetra = 0;
  pbf->quadrature_tetra = NULL;

  /* Set function pointers */
  if (dim == 3) {

    switch (k) {

    case 0:
      pbf->eval_all_at_point = _k0_eval_all_at_point;
      pbf->eval_at_point = _k0_eval_at_point;
      pbf->compute_projector = _ck0_compute_projector;
      pbf->compute_factorization = _k0_compute_facto;
      pbf->project = _k0_project;
      pbf->dump_projector = _cka_dump_projector;
      pbf->n_gpts_tetra = 4;
      pbf->quadrature_tetra = cs_quadrature_tet_4pts;
      break;
    case 1:
      pbf->eval_all_at_point = _ck1_eval_all_at_point;
      pbf->eval_at_point = _ck1_eval_at_point;
      pbf->compute_projector = _ck1_compute_projector;
      pbf->compute_factorization = _ck1_compute_facto;
      pbf->project = _ck1_project;
      pbf->dump_projector = _cka_dump_projector;
      pbf->n_gpts_tetra = 5;
      pbf->quadrature_tetra = cs_quadrature_tet_5pts;
      break;
    case 2:
      pbf->eval_all_at_point = _cka_eval_all_at_point;
      pbf->eval_at_point = _cka_eval_at_point;
      pbf->compute_projector = _ck2_compute_projector;
      pbf->compute_factorization = _ka_compute_facto;
      pbf->project = _ka_project;
      pbf->dump_projector = _cka_dump_projector;
      pbf->n_gpts_tetra = 15;
      pbf->quadrature_tetra = cs_quadrature_tet_15pts;
      break;
    default: /* Arbitrary order */
      pbf->eval_all_at_point = _cka_eval_all_at_point;
      pbf->eval_at_point = _cka_eval_at_point;
      pbf->compute_projector = _cka_compute_projector;
      pbf->compute_factorization = _ka_compute_facto;
      pbf->project = _ka_project;
      pbf->dump_projector = _cka_dump_projector;
      pbf->n_gpts_tetra = 15;
      pbf->quadrature_tetra = cs_quadrature_tet_15pts;
      break;

    }

    if (flag & CS_BASIS_FUNC_MONOMIAL)
      pbf->setup = _mono_cell_basis_setup;
    else  /* Rescaling based on the inertial axis */
      pbf->setup = _iner_cell_basis_setup;

  }
  else {

    switch (k) {

    case 0:
      pbf->eval_all_at_point = _k0_eval_all_at_point;
      pbf->eval_at_point = _k0_eval_at_point;
      pbf->compute_projector = _fk0_compute_projector;
      pbf->compute_factorization = _k0_compute_facto;
      pbf->project = _k0_project;
      pbf->dump_projector = _cka_dump_projector;
      pbf->n_gpts_tria = 3;
      pbf->quadrature_tria = cs_quadrature_tria_3pts;
      break;
    case 1:
      pbf->eval_all_at_point = _fk1_eval_all_at_point;
      pbf->eval_at_point = _fk1_eval_at_point;
      pbf->compute_projector = _fk1_compute_projector;
      pbf->compute_factorization = _fk1_compute_facto;
      pbf->project = _fk1_project;
      pbf->dump_projector = _fka_dump_projector;
      pbf->n_gpts_tria = 4;
      pbf->quadrature_tria = cs_quadrature_tria_4pts;
      break;
    case 2:
      pbf->eval_all_at_point = _fka_eval_all_at_point;
      pbf->eval_at_point = _fka_eval_at_point;
      pbf->compute_projector = _fk2_compute_projector;
      pbf->compute_factorization = _fk2_compute_facto;
      pbf->project = _fk2_project;
      pbf->dump_projector = _fka_dump_projector;
      pbf->n_gpts_tria = 7;
      pbf->quadrature_tria = cs_quadrature_tria_7pts;
      break;
    default: /* Arbitrary order */
      pbf->eval_all_at_point = _fka_eval_all_at_point;
      pbf->eval_at_point = _fka_eval_at_point;
      pbf->compute_projector = _fka_compute_projector;
      pbf->compute_factorization = _ka_compute_facto;
      pbf->project = _ka_project;
      pbf->dump_projector = _fka_dump_projector;
      pbf->n_gpts_tria = 7;
      pbf->quadrature_tria = cs_quadrature_tria_7pts;
      break;
    }

    if (flag & CS_BASIS_FUNC_MONOMIAL)
      pbf->setup = _mono_face_basis_setup;
    else
      pbf->setup = _iner_face_basis_setup;

  }

  return pbf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_basis_func_t structure which is associated to an
 *         existing set of basis functions.
 *         Up to now, only cell basis functions are handled.
 *         Building a projection matrix is not possible in this case.
 *
 * \param[in]  ref  set of basis function used as a reference
 *
 * \return a pointer to the new cs_basis_func_t for gradient of the current
 *         basis functions
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_grad_create(const cs_basis_func_t   *ref)
{
  assert(ref->dim == 3);

  /* Structure to handle polynomial basis functions */
  cs_basis_func_t  *gbf = NULL;

  BFT_MALLOC(gbf, 1, cs_basis_func_t);

  gbf->flag = ref->flag | CS_BASIS_FUNC_GRADIENT;
  gbf->poly_order = ref->poly_order; /* Grad(P^(k+1)_d) is of order k */
  gbf->dim = ref->dim;
  gbf->size = cs_math_binom(gbf->poly_order + 1 + ref->dim, ref->dim);
  gbf->phi0 = 1;  /* not useful */

  /* Copy axis and center */
  BFT_MALLOC(gbf->axis, ref->dim, cs_nvec3_t);

  /* Build a basis of polynomial functions of order k+1.
     Gradient is apply on-the-fly when the evaluation is performed */
  gbf->n_deg_elts = 0;
  gbf->deg = NULL;
  if (gbf->poly_order > 0) { /* Basis of order striclty greater than 1 */

    /* Remove the 4 first trivial basis functions in 3D */
    gbf->n_deg_elts = gbf->size -(ref->dim + 1);
    BFT_MALLOC(gbf->deg, gbf->n_deg_elts*ref->dim, short int);

    short int count = 0;
    for (short int ik = 2; ik < gbf->poly_order+2; ik++) {
      for (short int e1 = ik; e1 > -1; e1--) {
        for (short int e2 = ik-e1; e2 > -1; e2--) {
          gbf->deg[count*ref->dim    ] = e1;
          gbf->deg[count*ref->dim + 1] = e2;
          gbf->deg[count*ref->dim + 2] = ik - e1 -e2;
          count++;
        }
      }
    }
    assert(count == gbf->n_deg_elts);

  }

  /* Only the function pointer "eval_all_at_point" makes sense in the
   context of the gradient of a set of basis functions */
  gbf->setup = NULL; /* Copy from the reference. No need to set */
  gbf->projector = NULL;
  gbf->dump_projector = NULL;
  gbf->compute_projector = NULL;
  gbf->compute_factorization = NULL;
  gbf->facto_max_size = 0;
  gbf->facto = NULL;  /* array storing the matrix factorization */
  gbf->project = NULL;
  gbf->n_gpts_tria = ref->n_gpts_tria;
  gbf->quadrature_tria = ref->quadrature_tria;
  gbf->n_gpts_tetra = ref->n_gpts_tetra;
  gbf->quadrature_tetra = ref->quadrature_tetra;

  /* Set function pointers */
  switch (gbf->poly_order) {

  case 0:
    gbf->eval_all_at_point = _cgk1_eval_all_at_point;
    gbf->eval_at_point = _cgk1_eval_at_point;
    break;
  default: /* Arbitrary order */
    gbf->eval_all_at_point = _cgka_eval_all_at_point;
    gbf->eval_at_point = _cgka_eval_at_point;
    break;
  }

  return gbf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the center and the different axis from the reference basis
 *         Up to now, only cell basis functions are handled.
 *
 * \param[in]      ref   set of basis function used as a reference
 * \param[in, out] rcv   set of basis function where members are set
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_copy_setup(const cs_basis_func_t   *ref,
                         cs_basis_func_t         *rcv)
{
  assert(ref != NULL);

  /* Copy axis and center */
  for (int i = 0; i < ref->dim; i++) {
    for (int k = 0; k < 3; k++)
      rcv->axis[i].unitv[k] = ref->axis[i].unitv[k];
    rcv->axis[i].meas = ref->axis[i].meas;
  }
  for (int k = 0; k < 3; k++)
    rcv->center[k] = ref->center[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_basis_func_t structure
 *
 * \param[in, out]  pbf   pointer to the cs_basis_func_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_basis_func_t *
cs_basis_func_free(cs_basis_func_t  *pbf)
{
  if (pbf == NULL)
    return pbf;

  BFT_FREE(pbf->axis);
  BFT_FREE(pbf->deg);

  if (pbf->projector != NULL)
    pbf->projector = cs_sdm_free(pbf->projector);
  pbf->facto_max_size = 0;
  BFT_FREE(pbf->facto);

  BFT_FREE(pbf);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set options for basis functions when using HHO schemes
 *
 * \param[in]  face_flag    options related to face basis functinos
 * \param[in]  cell_flag    options related to cell basis functinos
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_set_hho_flag(cs_flag_t   face_flag,
                           cs_flag_t   cell_flag)
{
  cs_basis_func_hho_face_flag = face_flag;
  cs_basis_func_hho_cell_flag = cell_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get options for basis functions when using HHO schemes
 *
 * \param[out] face_flag   pointer to options related to face basis functinos
 * \param[out] cell_flag   pointer to options related to cell basis functinos
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_get_hho_flag(cs_flag_t   *face_flag,
                           cs_flag_t   *cell_flag)
{
  *face_flag = cs_basis_func_hho_face_flag;
  *cell_flag = cs_basis_func_hho_cell_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_basis_func_t structure
 *
 * \param[in]  pbf   pointer to the cs_basis_func_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_dump(const cs_basis_func_t  *pbf)
{
  cs_log_printf(CS_LOG_DEFAULT, "\n basis function: %p\n", (const void *)pbf);
  if (pbf == NULL)
    return;

  cs_log_printf(CS_LOG_DEFAULT,
                " flag: %d; dim; %d; poly_order: %d; size: %d\n",
                pbf->flag, pbf->dim, pbf->poly_order, pbf->size);
  cs_log_printf(CS_LOG_DEFAULT,
                " phi0: % .4e; center: (% .4e, % .4e % .4e)\n",
                pbf->phi0, pbf->center[0], pbf->center[1], pbf->center[2]);

  for (int k = 0; k < pbf->dim; k++)
    cs_log_printf(CS_LOG_DEFAULT,
                  " axis(%d) [% .4e, % .4e % .4e] % .4e\n",
                  k, pbf->axis[k].unitv[0], pbf->axis[k].unitv[1],
                  pbf->axis[k].unitv[2], pbf->axis[k].meas);

  if (pbf->deg != NULL) {
    for (int k = 0; k < pbf->dim; k++) {
      for (int i = 0; i < pbf->n_deg_elts; i++)
        cs_log_printf(CS_LOG_DEFAULT, "%3d", pbf->deg[i*pbf->dim+k]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a cs_basis_func_t structure
 *         Print into the file f if given otherwise open a new file named
 *         fname if given otherwise print into the standard output
 *
 * \param[in]  fp      pointer to a file structure or NULL
 * \param[in]  fname   filename or NULL
 * \param[in]  pbf     pointer to the cs_basis_func_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_basis_func_fprintf(FILE                   *fp,
                      const char             *fname,
                      const cs_basis_func_t  *pbf)
{
  FILE  *fout = stdout;
  if (fp != NULL)
    fout = fp;
  else if (fname != NULL) {
    fout = fopen(fname, "w");
  }

  fprintf(fout, "\n basis function: %p\n", (const void *)pbf);
  if (pbf == NULL)
    return;

  fprintf(fout, " flag: %d; dim; %d; poly_order: %d; size: %d\n",
          pbf->flag, pbf->dim, pbf->poly_order, pbf->size);
  fprintf(fout, " phi0: % .4e; center: (% .4e, % .4e % .4e)\n",
          pbf->phi0, pbf->center[0], pbf->center[1], pbf->center[2]);

  for (int k = 0; k < pbf->dim; k++)
    fprintf(fout, " axis(%d) [% .5e, % .5e % .5e] % .4e\n",
            k, pbf->axis[k].unitv[0], pbf->axis[k].unitv[1],
            pbf->axis[k].unitv[2], pbf->axis[k].meas);

  if (pbf->deg != NULL) {
    for (int k = 0; k < pbf->dim; k++) {
      for (int i = 0; i < pbf->n_deg_elts; i++)
        fprintf(fout, "%3d", pbf->deg[i*pbf->dim+k]);
      fprintf(fout, "\n");
    }
  }

  if (pbf->facto != NULL) {
    const int facto_size = ((pbf->size + 1)*pbf->size)/2;
    fprintf(fout, "Factorization:\n");
    for (int i = 0; i < facto_size; i++)
      fprintf(fout, " % -9.5e", pbf->facto[i]);
    fprintf(fout, "\n");
  }

  if (fout != stdout && fout != fp)
    fclose(fout);
}

/*----------------------------------------------------------------------------*/

#undef _dp3
#undef CK1_SIZE
#undef CK2_SIZE
#undef FK1_SIZE
#undef FK2_SIZE

END_C_DECLS
