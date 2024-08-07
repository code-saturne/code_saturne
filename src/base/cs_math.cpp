/*============================================================================
 * Mathematical base functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_math.c
           Mathematical base functions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Numerical constants */

#if !(defined(__NVCC__) && defined(__CUDA_ARCH__))

const cs_real_t cs_math_zero_threshold = FLT_MIN;
const cs_real_t cs_math_1ov3 = 1./3.;
const cs_real_t cs_math_2ov3 = 2./3.;
const cs_real_t cs_math_4ov3 = 4./3.;
const cs_real_t cs_math_5ov3 = 5./3.;
const cs_real_t cs_math_1ov6 = 1./6.;
const cs_real_t cs_math_1ov12 = 1./12.;
const cs_real_t cs_math_1ov24 = 1./24.;

/*! epsilon \f$ 10^{-12}\f$ */
const cs_real_t cs_math_epzero = 1e-12;

/*! infinite \f$ 10^{+30}\f$ */
const cs_real_t cs_math_infinite_r = 1.e30;

/*! big value \f$ 10^{+12}\f$ */
const cs_real_t cs_math_big_r = 1.e12;

/*! \f$ \pi \f$ value with 20 digits */
const cs_real_t cs_math_pi = 3.14159265358979323846;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a rotation of the elements in a 3x3 real values matrix
 *
 * \param[in, out]  m   matrix of 3x3 real values
 * \param[in]       i   1st index
 * \param[in]       j   2nd index
 * \param[in]       k   3rd index
 * \param[in]       l   4th index
 * \param[in]       s   rate
 * \param[in]       t   rate
 */
/*----------------------------------------------------------------------------*/

static inline void
_rotate_ind_33(cs_real_t  m[3][3],
               cs_lnum_t  i,
               cs_lnum_t  j,
               cs_lnum_t  k,
               cs_lnum_t  l,
               cs_real_t  s,
               cs_real_t  t)
{
  /* Save values of m[i][j] and m[k][l] */
  cs_real_t m_ij = m[i][j];
  cs_real_t m_kl = m[k][l];

  /* Modify the values of (i,j) and (k,l) */
  m[i][j] = m_ij - s*(m_kl + m_ij*t);
  m[k][l] = m_kl + s*(m_ij - m_kl*t);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_math_sym_33_inv_cramer(const cs_real_t s[6],
                            cs_real_t       sout[6]);

void
cs_f_math_sym_33_product(const cs_real_t  s1[6],
                         const cs_real_t  s2[6],
                         cs_real_t        sout[6]);

void
cs_f_math_reduce_sym_prod_33_to_66(const cs_real_t       s[3][3],
                                   cs_real_t  (*restrict sout)[6]);

void
cs_f_math_3_normalize(const cs_real_t vin[3],
                      cs_real_t       vout[3]);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_inv_cramer
 *----------------------------------------------------------------------------*/

void
cs_f_math_sym_33_inv_cramer(const cs_real_t s[6],
                            cs_real_t       sout[6])
{
  cs_math_sym_33_inv_cramer(s, sout);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_product
 *----------------------------------------------------------------------------*/

void
cs_f_math_sym_33_product(const cs_real_t  s1[6],
                         const cs_real_t  s2[6],
                         cs_real_t        sout[6])
{
  cs_math_sym_33_product(s1, s2, sout);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_reduce_sym_prod_33_to_66
 *----------------------------------------------------------------------------*/

void
cs_f_math_reduce_sym_prod_33_to_66(const cs_real_t       s[3][3],
                                   cs_real_t  (*restrict sout)[6])
{
  cs_math_reduce_sym_prod_33_to_66(s, sout);
}

/*----------------------------------------------------------------------------
 * Wrapper
 *----------------------------------------------------------------------------*/

void
cs_f_math_3_normalize(const cs_real_t vin[3],
                      cs_real_t       vout[3])
{
  cs_math_3_normalize(vin, vout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute all eigenvalues of a 3x3 symmetric matrix
 *        with symmetric storage.
 *
 * Based on: Oliver K. Smith "eigenvalues of a symmetric 3x3 matrix",
 *           Communication of the ACM (April 1961)
 *           (Wikipedia article entitled "Eigenvalue algorithm")
 *
 * \param[in]  m          3x3 symmetric matrix (m11, m22, m33, m12, m23, m13)
 * \param[out] eig_vals   size 3 vector
 */
/*----------------------------------------------------------------------------*/

void
cs_math_sym_33_eigen(const cs_real_t  m[6],
                     cs_real_t        eig_vals[3])
{
  constexpr cs_real_t c_1ov3 = 1./3.;
  constexpr cs_real_t c_1ov6 = 1./6.;

  cs_real_t  e, e1, e2, e3;

  cs_real_t  p1 = cs_math_3_square_norm((const cs_real_t *)(m+3));
  cs_real_t  d2 = cs_math_3_square_norm((const cs_real_t *)m);

  cs_real_t  tr = (m[0] + m[1] + m[2]);
  cs_real_t  tr_third = c_1ov3 * tr;

  e1 = m[0] - tr_third, e2 = m[1] - tr_third, e3 = m[2] - tr_third;
  cs_real_t  p2 = e1*e1 + e2*e2 + e3*e3 + 2.*p1;

  cs_real_t  p = sqrt(p2*c_1ov6);

  if (p1 > cs_math_epzero*d2 && p > 0.) { /* m is not diagonal */

    cs_real_6_t  n;
    cs_real_t  ovp = 1./p;

    for (int  i = 0; i < 3; i++) {
      /* Diagonal */
      n[i] = ovp * (m[i] - tr_third);
      /* Extra diagonal */
      n[3 + i] = ovp * m[3 + i];
    }

    /* r should be between -1 and 1 but truncation error and bad conditioning
       can lead to slightly under/over-shoot */
    cs_real_t  r = 0.5 * cs_math_sym_33_determinant(n);

    cs_real_t  cos_theta, cos_theta_2pi3;
    if (r <= -1.) {
      cos_theta = 0.5; /* theta = pi/3; */
      cos_theta_2pi3 = -1.;
    }
    else if (r >= 1.) {
      cos_theta = 1.; /* theta = 0.; */
      cos_theta_2pi3 = -0.5;
    }
    else {
      cos_theta = cos(c_1ov3*acos(r));
      cos_theta_2pi3 = cos(c_1ov3*(acos(r) + 2.*cs_math_pi));
    }

    /* eigenvalues computed should satisfy e1 < e2 < e3 */
    e3 = tr_third + 2.*p*cos_theta;
    e1 = tr_third + 2.*p*cos_theta_2pi3;
    e2 = tr - e1 -e3; // since tr(m) = e1 + e2 + e3

  }
  else { /* m is diagonal */

    e1 = m[0], e2 = m[1], e3 = m[2];

  } /* diagonal or not */

  if (e3 < e2) e = e3, e3 = e2, e2 = e;
  if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
  else {
    if (e2 < e1) e = e2, e2 = e1, e1 = e;
  }

  /* Return values */
  eig_vals[0] = e1;
  eig_vals[1] = e2;
  eig_vals[2] = e3;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute max/min eigenvalues ratio and max. eigenvalue of a 3x3
 *         symmetric matrix with non-symmetric storage.
 *
 * Based on: Oliver K. Smith "eigenvalues of a symmetric 3x3 matrix",
 *           Communication of the ACM (April 1961)
 *           (Wikipedia article entitled "Eigenvalue algorithm")
 *
 * \param[in]  m          3x3 matrix
 * \param[out] eig_ratio  max/min
 * \param[out] eig_max    max. eigenvalue
 */
/*----------------------------------------------------------------------------*/

void
cs_math_33_eigen(const cs_real_t     m[3][3],
                 cs_real_t          *eig_ratio,
                 cs_real_t          *eig_max)
{
  constexpr cs_real_t c_1ov3 = 1./3.;
  constexpr cs_real_t c_1ov6 = 1./6.;

  cs_real_t  e, e1, e2, e3;

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
  e1 = m[0][1]-m[1][0], e2 = m[0][2]-m[2][0], e3 = m[1][2]-m[2][1];
  if (e1*e1 + e2*e2 + e3*e3 > 0.0)
    bft_error(__FILE__, __LINE__, 0,
              " The given 3x3 matrix is not symmetric.\n"
              " Stop computing eigenvalues of a 3x3 matrix since the"
              " algorithm is only dedicated to symmetric matrix.");
#endif

  cs_real_t  p1 = m[0][1]*m[0][1] + m[0][2]*m[0][2] + m[1][2]*m[1][2];

  if (p1 > 0.0) { /* m is not diagonal */

    cs_real_t  theta;
    cs_real_t  n[3][3];
    cs_real_t  tr = c_1ov3*(m[0][0] + m[1][1] + m[2][2]);

    e1 = m[0][0] - tr, e2 = m[1][1] - tr, e3 = m[2][2] - tr;
    cs_real_t  p2 = e1*e1 + e2*e2 + e3*e3 + 2*p1;

    assert(p2 > 0);
    cs_real_t  p = sqrt(p2*c_1ov6);
    cs_real_t  ovp = 1./p;

    for (int  i = 0; i < 3; i++) {
      n[i][i] = ovp * (m[i][i] - tr);
      for (int j = i + 1; j < 3; j++) {
        n[i][j] = ovp*m[i][j];
        n[j][i] = n[i][j];
      }
    }

    /* r should be between -1 and 1 but truncation error and bad conditioning
       can lead to slightly under/over-shoot */
    cs_real_t  r = 0.5 * cs_math_33_determinant((const cs_real_t (*)[3])n);

    if (r <= -1)
      theta = c_1ov3*cs_math_pi;
    else if (r >= 1)
      theta = 0.;
    else
      theta = c_1ov3*acos(r);

    /* eigenvalues computed should satisfy e1 < e2 < e3 */
    e3 = tr + 2*p*cos(theta);
    e1 = tr + 2*p*cos(theta + 2*cs_math_pi*c_1ov3);
    e2 = 3*tr - e1 -e3; // since tr(m) = e1 + e2 + e3

  }
  else { /* m is diagonal */

    e1 = m[0][0], e2 = m[1][1], e3 = m[2][2];
    if (e3 < e2) e = e3, e3 = e2, e2 = e;
    if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
    else {
      if (e2 < e1) e = e2, e2 = e1, e1 = e;
    }

  } /* diagonal or not */

  /* Return values */
  if (fabs(e1) > 0)
    *eig_ratio = e3/e1;
  else
    *eig_ratio = 1;
  *eig_max = e3;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the length (Euclidean norm) between two points xa and xb in
 *         a Cartesian coordinate system of dimension 3
 *
 * \param[in]   xa       coordinate of the first extremity
 * \param[in]   xb       coordinate of the second extremity
 * \param[out]  len      pointer to the length of the vector va -> vb
 * \param[out]  unitv    unitary vector along va -> vb
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE void
cs_math_3_length_unitv(const cs_real_t    xa[3],
                       const cs_real_t    xb[3],
                       cs_real_t         *len,
                       cs_real_3_t        unitv)
{
  cs_real_t  invl;
  cs_real_3_t  diff;

  diff[0] = xb[0] - xa[0];
  diff[1] = xb[1] - xa[1];
  diff[2] = xb[2] - xa[2];

  *len = cs_math_3_norm(diff), invl = 1/(*len);

  unitv[0] = invl*diff[0];
  unitv[1] = invl*diff[1];
  unitv[2] = invl*diff[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the convex_hull generated by 3 points.
 *         This corresponds to the computation of the surface of a triangle
 *
 * \param[in]  xv  coordinates of the first vertex
 * \param[in]  xe  coordinates of the second vertex
 * \param[in]  xf  coordinates of the third vertex
 *
 * \return the surface of a triangle
 */
/*----------------------------------------------------------------------------*/

double
cs_math_surftri(const cs_real_t  xv[3],
                const cs_real_t  xe[3],
                const cs_real_t  xf[3])
{
  cs_real_t u[3], v[3], cp[3];

  for (int k = 0; k < 3; k++) {
    u[k] = xe[k] - xv[k];
    v[k] = xf[k] - xv[k];
  }
  cs_math_3_cross_product(u, v, cp);

  return  0.5 * cs_math_3_norm(cp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume of the convex_hull generated by 4 points.
 *         This is equivalent to the computation of the volume of a tetrahedron
 *
 * \param[in]  xv  coordinates of the first vertex
 * \param[in]  xe  coordinates of the second vertex
 * \param[in]  xf  coordinates of the third vertex
 * \param[in]  xc  coordinates of the fourth vertex
 *
 * \return the volume of the tetrahedron.
 */
/*----------------------------------------------------------------------------*/

double
cs_math_voltet(const cs_real_t   xv[3],
               const cs_real_t   xe[3],
               const cs_real_t   xf[3],
               const cs_real_t   xc[3])
{
  constexpr cs_real_t c_1ov6 = 1./6.;

  double  lev, lef, lec;
  cs_real_3_t  uev, uef, uec, ucp;

  cs_math_3_length_unitv(xe, xv, &lev, uev);
  cs_math_3_length_unitv(xe, xf, &lef, uef);
  cs_math_3_length_unitv(xe, xc, &lec, uec);
  cs_math_3_cross_product(uev, uef, ucp);

  return  c_1ov6 *lev*lef*lec* fabs(cs_math_3_dot_product(ucp, uec));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate eigenvalues and eigenvectors
 *         of a real symmetric matrix m1[3,3]: m1*m2 = lambda*m2
 *
 * Use of Jacobi method for symmetric matrices
 * (adapted from the book Numerical Recipes in C, Chapter 11.1)
 *
 * \param[in]     m_in         matrix of 3x3 real values (initial)
 * \param[in]     tol_err      absolute tolerance (sum of off-diagonal elements)
 * \param[out]    eig_val      vector of 3 real values (eigenvalues)
 * \param[out]    eig_vec      matrix of 3x3 real values (eigenvectors)
 */
/*----------------------------------------------------------------------------*/

void
cs_math_33_eig_val_vec(const cs_real_t   m_in[3][3],
                       const cs_real_t   tol_err,
                       cs_real_t        *restrict eig_val,
                       cs_real_t       (*restrict eig_vec)[3])
{
  /* Declaration of local variables
   * vec1, vec2:  vectors of 3 real values (copies of diagonal values)
   * m:           matrix of 3x3 real values (copy of m_in)
   * epsilon:     error (like machine epsilon for floats) */
  cs_real_t vec1[3], vec2[3];
  cs_real_t m[3][3] = {{m_in[0][0], m_in[0][1], m_in[0][2]},
                       {m_in[1][0], m_in[1][1], m_in[1][2]},
                       {m_in[2][0], m_in[2][1], m_in[2][2]}};
  cs_real_t epsilon = 1.0e-16;

  for (int id = 0; id < 3; id++) {
    eig_val[id] = m_in[id][id];
    vec1[id] = eig_val[id];
    vec2[id] = 0.0;
  }

  /* The strategy here is to adopt a cyclic Jacobi method
   * where the diagonalisation is carried out by several sweeps
   * (each sweep to increase the precision of the result)
   * Here, we perform up to 50 sweeps (the algorithm stops when
   * reaching the given tolerance) */

  // Error = sum off-diagonal elements
  cs_real_t error
    = cs_math_fabs(m[0][1]) + cs_math_fabs(m[0][2]) + cs_math_fabs(m[1][2]);

  // Iterative solution
  for (int i_sweep = 0; (i_sweep < 50) && (error > tol_err); i_sweep++) {

    // Start loop on off-diagonal elements
    for (int id1 = 0; id1 < 2; id1 ++) {
      for (int id2 = id1+1; id2 < 3; id2 ++) {
        // After 4 sweeps, skip rotation if off-diagonal element is small
        if (cs_math_fabs(m[id1][id2]) < epsilon) {
          if (i_sweep > 4)
            m[id1][id2] = 0.0;
        }
        // Otherwise, ...
        else { // if (abs(m[id1][id2]) >= epsilon)
          cs_real_t val1, val2, val3, val4, val5;    // Declare five variables val.
          // Get val1
          cs_real_t diff = eig_val[id2] - eig_val[id1];
          if (cs_math_fabs(m[id1][id2]) < epsilon)
            val1 = m[id1][id2] / diff;
          else {
            cs_real_t theta = 0.5 * diff / m[id1][id2];
            val1 = 1.0 / (cs_math_fabs(theta)+sqrt(1.0+cs_math_pow2(theta)));
            if ( theta < 0 )
              val1 = -val1;
          }
          // Get val2, val3 and val4
          val3 = 1.0 / sqrt(1.0+cs_math_pow2(val1));
          val2 = val1*val3;
          val4 = val2 / (1.0 + val3);
          val5 = val1 * m[id1][id2];
          // Accumulate correction to diagonal elements
          vec2[id1] -= val5;
          vec2[id2] += val5;
          eig_val[id1] -= val5;
          eig_val[id2] += val5;

          m[id1][id2] = 0.0;
          // Rotate
          for (int id3 = 0; id3 <= id1-1; id3++) // Rotations 0 <= j < p
            _rotate_ind_33(m, id3, id1, id3, id2, val2, val4);
          for (int id3 = id1+1; id3 <= id2-1; id3++) // Rotations p < j < q
            _rotate_ind_33(m, id1, id3, id3, id2, val2, val4);
          for (int id3 = id2+1; id3 < 3; id3++) // Rotations q < j <= n
            _rotate_ind_33(m, id1, id3, id2, id3, val2, val4);
          for (int id3 = 0; id3 < 3; id3++)
            _rotate_ind_33(eig_vec, id3, id1, id3, id2, val2, val4);
        }
      }
    }

    // Update d and reinitialize z
    for (int id = 0; id < 3; id++ ) {
      vec1[id] += vec2[id];
      eig_val[id] = vec1[id];
      vec2[id] = 0.0;
    }

    // Update the error
    error = cs_math_fabs(m[0][1]) + cs_math_fabs(m[0][2]) + cs_math_fabs(m[1][2]);
  }

  /* Sort eigenvalues and eigenvectors with ascending order */
  for (int id1 = 0; id1 < 2; id1++) {
    cs_lnum_t ind_min = id1;
    for (int id2 = id1+1; id2 < 3; id2++) {
      if ( eig_val[id2] < eig_val[ind_min] )
        ind_min = id2;
    }
    if ( ind_min != id1 ) {
      cs_real_t temp = eig_val[ind_min];
      eig_val[ind_min] = eig_val[id1];
      eig_val[id1] = temp;
      for (int id2 = 0; id2 < 3; id2++) {
        temp = eig_vec[id2][ind_min];
        eig_vec[id2][ind_min] = eig_vec[id2][id1];
        eig_vec[id2][id1] = temp;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute LU factorization of an array of dense matrices
 *        of identical size.
 *
 * \param[in]   n_blocks  number of blocks
 * \param[in]   b_size    block size
 * \param[in]   a         matrix blocks
 * \param[out]  a_lu      LU factorizations of matrix blocks
 */
/*----------------------------------------------------------------------------*/

void
cs_math_fact_lu(cs_lnum_t         n_blocks,
                int               b_size,
                const cs_real_t  *a,
                cs_real_t        *a_lu)
{
# pragma omp parallel for if(n_blocks > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_blocks; i++) {

    cs_real_t *restrict _a_lu = &a_lu[b_size*b_size*i];
    const cs_real_t *restrict  _a = &a[b_size*b_size*i];

    _a_lu[0] = _a[0];
    for (cs_lnum_t ii = 1; ii < b_size; ii++) {
      _a_lu[ii] = _a[ii];
      _a_lu[ii*b_size] = _a[ii*b_size]/_a[0];
    }
    for (cs_lnum_t ii = 1; ii < b_size - 1; ii++) {
      _a_lu[ii + ii*b_size] = _a[ii + ii*b_size];
      for (cs_lnum_t kk = 0; kk < ii; kk++) {
        _a_lu[ii + ii*b_size] -= _a_lu[ii*b_size + kk]
                                *_a_lu[kk*b_size + ii];
      }

      for (cs_lnum_t jj = ii + 1; jj < b_size; jj++) {
        _a_lu[ii*b_size + jj] = _a[ii*b_size + jj];
        _a_lu[jj*b_size + ii] =   _a[jj*b_size + ii]
                                / _a_lu[ii*b_size + ii];
        for (cs_lnum_t kk = 0; kk < ii; kk++) {
          _a_lu[ii*b_size + jj] -=  _a_lu[ii*b_size + kk]
                                   *_a_lu[kk*b_size + jj];
          _a_lu[jj*b_size + ii] -=  _a_lu[jj*b_size + kk]
                                   *_a_lu[kk*b_size + ii]
                                   /_a_lu[ii*b_size + ii];
        }
      }
    }
    _a_lu[b_size*b_size -1] = _a[b_size*b_size - 1];

    for (cs_lnum_t kk = 0; kk < b_size - 1; kk++) {
      _a_lu[b_size*b_size - 1] -=  _a_lu[(b_size-1)*b_size + kk]
                                  *_a_lu[kk*b_size + b_size -1];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Block Jacobi utilities.
 *         Compute forward and backward to solve an LU P*P system.
 *
 * \param[in]   a_lu   matrix LU factorization
 * \param[in]   n      matrix size
 * \param[out]  x      solution
 * \param[out]  b      right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_math_fw_and_bw_lu(const cs_real_t  a_lu[],
                     int              n,
                     cs_real_t        x[],
                     const cs_real_t  b[])
{
  cs_real_t  _aux[256];
  cs_real_t  *aux = _aux;

  if (n > 256)
    BFT_MALLOC(aux, n, cs_real_t);

  /* forward */
  for (int ii = 0; ii < n; ii++) {
    aux[ii] = b[ii];
    for (int jj = 0; jj < ii; jj++) {
      aux[ii] -= aux[jj]*a_lu[ii*n + jj];
    }
  }

  /* backward */
  for (int ii = n - 1; ii >= 0; ii-=1) {
    x[ii] = aux[ii];
    for (int jj = n - 1; jj > ii; jj-=1) {
      x[ii] -= x[jj]*a_lu[ii*n + jj];
    }
    x[ii] /= a_lu[ii*(n + 1)];
  }

  if (n > 256)
    BFT_FREE(aux);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (m00, m10, m11, m20, m21, m22, m30, m31, m32, m33) in
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33) out
 */
/*----------------------------------------------------------------------------*/

void
cs_math_sym_44_factor_ldlt(cs_real_t  ldlt[10])
{
  /* Factorization */

  // j=0
  const cs_real_t  d00 = ldlt[0]; // m00
  assert(fabs(d00) > cs_math_zero_threshold);

  const cs_real_t  f00 = 1. / d00;
  const cs_real_t  l10 = ldlt[1] * f00; // m01
  const cs_real_t  l20 = ldlt[3] * f00; // m02
  const cs_real_t  l30 = ldlt[6] * f00; // m03

  // j=1
  const cs_real_t  d11 = ldlt[2] - l10*l10*d00; // m11
  assert(fabs(d11) > cs_math_zero_threshold);
  const cs_real_t  f11 = 1. / d11;
  const cs_real_t  l21 = (ldlt[4] - l20*d00*l10) * f11; // m12
  const cs_real_t  l31 = (ldlt[7] - l30*d00*l10) * f11; // m13

  // j=2
  const cs_real_t  d22 = ldlt[5] - l20*d00*l20 - l21*d11*l21;  // m22
  assert(fabs(d22) > cs_math_zero_threshold);
  const cs_real_t  f22 = 1. / d22;
  const cs_real_t  l32 = (ldlt[8] - l30*d00*l20 - l31*d11*l21) * f22;  // m23

  // j=3
  const cs_real_t  d33 = ldlt[9] - l30*d00*l30 - l31*d11*l31 - l32*d22*l32; // m33
  assert(fabs(d33) > cs_math_zero_threshold);
  const cs_real_t  f33 = 1. / d33;

  ldlt[0] = f00;
  ldlt[1] = l10;
  ldlt[2] = f11;
  ldlt[3] = l20;
  ldlt[4] = l21;
  ldlt[5] = f22;
  ldlt[6] = l30;
  ldlt[7] = l31;
  ldlt[8] = l32;
  ldlt[9] = f33;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * Here we only need to use the last element of the solution vector,
 * so we return that value only and simplify the computation.
 *
 * \param[in, out]  ldlt  pointer to matrix coefficients:
 *                        (f00, l10, f11, l20, l21, f22, l30, l31, l32, f33)
 * \param[in]       rhs   right-hand side
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_math_sym_44_partial_solve_ldlt(const cs_real_t  ldlt[10],
                                  const cs_real_t  rhs[4])
{
  /* f00, f11, f22, f33, l32, l31, l30, l21, l20, l10
     0    1    2    3    4    5    6    7    8    9   */

  cs_real_t x[4]; /* solution */

  x[0] = rhs[0];
  x[1] = rhs[1] - x[0]*ldlt[1];
  x[2] = rhs[2] - x[0]*ldlt[3] - x[1]*ldlt[4];
  x[3] = rhs[3] - x[0]*ldlt[6] - x[1]*ldlt[7] - x[2]*ldlt[8];

  x[3] = x[3]*ldlt[9];

  return x[3];

  /*
    x[2] = x[2]*ldlt[5] - ldlt[8]*x[3];
    x[1] = x[1]*ldlt[2] - ldlt[7]*x[3] - ldlt[4]*x[2];
    x[0] = x[0]*ldlt[0] - ldlt[6]*x[3] - ldlt[3]*x[2] - ldlt[1]*x[1];
  */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
