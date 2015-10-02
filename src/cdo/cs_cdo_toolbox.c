/*============================================================================
 * Basic operations: dot product, cross product, sum...
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

struct _subsum_t {

  int      size;
  int     *idx;
  double  *sums;

};

static struct _subsum_t  _op_subsum;

static const double  _oversix = 1/6.0;
static const double  _overdim = 1/3.0;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Intialize by default a cs_data_info_t structure according to the
 *          datatype
 *
 * \param[in]  datatype
 *
 * \return a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_data_info_t
_init_dinfo(cs_datatype_t   datatype)
{
  cs_data_info_t  info;

  info.mean = 0.0;
  info.sigma = 0.0;
  info.euclidean_norm = 0.0;

  switch (datatype) {

  case CS_DOUBLE:
    info.min.value = DBL_MAX;
    info.max.value = -DBL_MAX;
    break;

  case CS_INT32:
    info.min.number = INT_MAX;
    info.max.number = -INT_MAX;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute simple information about an array of data.
 *          >> Algorithm from Mark Hoemmen (U.C. Berkeley)
 *
 * \param[in]      n_elts    number of couples in data
 * \param[in]      data      buffer containing input data
 * \param[in, out] info      pointer to a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_info_double(cs_lnum_t         n_elts,
                     const cs_real_t   data[],
                     cs_data_info_t   *info)
{
  int  i;

  if (n_elts == 0)
    return;

  /* Compute statistics */
  for (i = 0; i < n_elts; i++) {

    cs_real_t  val = data[i];
    cs_real_t  delta = val - info->mean;
    cs_real_t  chi = delta/(i+1);

    if (info->min.value > val)  info->min.value = val;
    if (info->max.value < val)  info->max.value = val;

    info->sigma += i*chi*delta;
    info->mean += chi;

  }
  info->sigma = sqrt(fabs(info->sigma)/n_elts);

  info->euclidean_norm = cs_euclidean_norm(n_elts, data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute simple information about an array of data.
 *          >> Algorithm from Mark Hoemmen (U.C. Berkeley)
 *
 * \param[in]      n_elts    number of couples in data
 * \param[in]      data      buffer containing input data
 * \param[in, out] info      pointer to a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_info_int32(cs_lnum_t         n_elts,
                    const cs_lnum_t   data[],
                    cs_data_info_t   *info)
{
  int  i;

  if (n_elts == 0)
    return;

  /* Compute statistics */
  for (i = 0; i < n_elts; i++) {

    cs_lnum_t  val = data[i];
    size_t  val2 = val*val;
    double  delta = val - info->mean;
    double  chi = delta/(i+1);

    if (info->min.number > val)  info->min.number = val;
    if (info->max.number < val)  info->max.number = val;

    info->sigma += i*chi*delta;
    info->mean += chi;
    info->euclidean_norm += val2;

  }

  info->sigma = sqrt(fabs(info->sigma)/n_elts);
  info->euclidean_norm = sqrt(fabs(info->euclidean_norm));
}

/*============================================================================
 * Inline public function prototypes for frequent usage
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute a dot product for vector of dimension 3
 *
 * \param[in]  u     first vector
 * \param[in]  v     second vector
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

inline double
_dp3(const cs_real_3_t  u,
     const cs_real_3_t  v)
{
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm of a vector of dimension 3
 *
 * \param[in]  v
 *
 * \return the norm value
 */
/*----------------------------------------------------------------------------*/

inline double
_n3(const cs_real_3_t  v)
{
  return sqrt(_dp3(v, v));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the length (euclidien norm) between two points xa and xb in
 *         a cartesian coordinate system of dimension 3
 *
 * \param[in]  xa   first coordinate
 * \param[in]  xb   second coordinate
 *
 * \return the length (in euclidean norm) between two points xa and xb
 */
/*----------------------------------------------------------------------------*/

inline double
_length3(const cs_real_3_t  xa,
         const cs_real_3_t  xb)
{
  cs_real_3_t  diff;

  diff[0] = xb[0] - xa[0];
  diff[1] = xb[1] - xa[1];
  diff[2] = xb[2] - xa[2];

  return _n3(diff);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the length (euclidien norm) between two points xa and xb in
 *         a cartesian coordinate system of dimension 3
 *
 * \param[in]   xa       coordinate of the first extremity
 * \param[in]   xb       coordinate of the second extremity
 * \param[out]  len      pointer to the length of the vector va -> vb
 * \param[out]  unit     unitary vector anlong va -> vb
 */
/*----------------------------------------------------------------------------*/

inline void
_lenunit3(const cs_real_3_t   xa,
          const cs_real_3_t   xb,
          cs_real_t          *len,
          cs_real_3_t        *unit)
{
  cs_real_t  invl;
  cs_real_3_t  diff;

  diff[0] = xb[0] - xa[0];
  diff[1] = xb[1] - xa[1];
  diff[2] = xb[2] - xa[2];
  *len = _n3(diff), invl = 1/(*len);
  unit[0][0] = invl*diff[0];
  unit[0][1] = invl*diff[1];
  unit[0][2] = invl*diff[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cross product of a vector of dimension 3
 *
 * \param[in]   u   first vector
 * \param[in]   v   second vector
 * \param[out]  w   result of u x v
 */
/*----------------------------------------------------------------------------*/

inline void
_cp3(const cs_real_3_t   u,
     const cs_real_3_t   v,
     cs_real_3_t        *w)
{
  w[0][0] = u[1]*v[2] - u[2]*v[1];
  w[0][1] = u[2]*v[0] - u[0]*v[2];
  w[0][2] = u[0]*v[1] - u[1]*v[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the 3x3 matrice by vector product
 *
 * \param[in]      m    a 3x3 matrix
 * \param[in]      v    a vector
 * \param[in, out] mv   pointer to the vector resulting of the matrix-vector op.
 */
/*----------------------------------------------------------------------------*/

inline void
_mv3(const cs_real_33_t  m,
     const cs_real_3_t   v,
     cs_real_3_t        *mv)
{
  mv[0][0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  mv[0][1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  mv[0][2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the determinant of a 3x3 matrix
 *
 * \param[in]  m    matrix
 *
 * \return the determinant
 */
/*----------------------------------------------------------------------------*/

inline cs_real_t
_detmat33(const cs_real_33_t   m)
{
  cs_real_t  com0 = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  cs_real_t  com1 = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  cs_real_t  com2 = m[0][1]*m[1][2] - m[1][1]*m[0][2];

  return m[0][0]*com0 + m[1][0]*com1 + m[2][0]*com2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 matrix
 *
 * \param[in]  in    matrix to inverse
 * \param[out] inv   inversed matrix
 */
/*----------------------------------------------------------------------------*/

inline void
_invmat33(const cs_real_33_t   in,
          cs_real_33_t        *inv)
{
  int  k, l;
  double  det, invdet;

  inv[0][0][0] = in[1][1]*in[2][2] - in[2][1]*in[1][2];
  inv[0][0][1] = in[2][1]*in[0][2] - in[0][1]*in[2][2];
  inv[0][0][2] = in[0][1]*in[1][2] - in[1][1]*in[0][2];

  inv[0][1][0] = in[2][0]*in[1][2] - in[1][0]*in[2][2];
  inv[0][1][1] = in[0][0]*in[2][2] - in[2][0]*in[0][2];
  inv[0][1][2] = in[1][0]*in[0][2] - in[0][0]*in[1][2];

  inv[0][2][0] = in[1][0]*in[2][1] - in[2][0]*in[1][1];
  inv[0][2][1] = in[2][0]*in[0][1] - in[0][0]*in[2][1];
  inv[0][2][2] = in[0][0]*in[1][1] - in[1][0]*in[0][1];

  det = in[0][0]*inv[0][0][0] + in[1][0]*inv[0][0][1] + in[2][0]*inv[0][0][2];
  assert(fabs(det) > DBL_MIN*1e3);  /* inversibility ? */
  invdet = 1 / det;

  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
      inv[0][k][l] *= invdet;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the eigenvalues of a 3x3 matrix which is symmetric and real
 *         -> Oliver K. Smith "eigenvalues of a symmetric 3x3 matrix",
 *         Communication of the ACM (April 1961)
 *         -> Wikipedia article entitled "Eigenvalue algorithm"
 *
 * \param[in]  m          3x3 matrix
 * \param[out] eig_ratio  max/min
 * \param[out] eig_max    max. eigenvalue
 */
/*----------------------------------------------------------------------------*/

void
cs_eigen_mat33(const cs_real_33_t  m,
               cs_real_t          *eig_ratio,
               cs_real_t          *eig_max)
{
  cs_real_t  e, e1, e2, e3;

  /* Sanity check */
  e1 = m[0][1]-m[1][0], e2 = m[0][2]-m[2][0], e3 = m[1][2]-m[2][1];
  assert(e1*e1 + e2*e2 + e3*e3 <= 0.0);

  cs_real_t  p1 = m[0][1]*m[0][1] + m[0][2]*m[0][2] + m[1][2]*m[1][2];

  if (p1 <= 0.0) { // m is diagonal

    e1 = m[0][0], e2 = m[1][1], e3 = m[2][2];
    if (e3 < e2) e = e3, e3 = e2, e2 = e;
    if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
    else {
      if (e2 < e1) e = e2, e2 = e1, e1 = e;
    }

  }
  else { // m is not diagonal

    cs_real_t  theta;
    cs_real_33_t  n;

    cs_real_t  tr = _overdim*(m[0][0] + m[1][1] + m[2][2]);

    e1 = m[0][0] - tr, e2 = m[1][1] - tr, e3 = m[2][2] - tr;
    cs_real_t  p2 = e1*e1 + e2*e2 + e3*e3 + 2*p1;

    assert(p2 > 0);
    cs_real_t  p = sqrt(p2*_oversix);
    cs_real_t  ovp = 1./p;

    for (int  i = 0; i < 3; i++) {
      n[i][i] = ovp * (m[i][i] - tr);
      for (int j = i + 1; j < 3; j++) {
        n[i][j] = ovp*m[i][j];
        n[j][i] = n[i][j];
      }
    }

    /* r should be between -1 and 1 but truncation error and bad conditionning
       can lead to slighty under/over-shoot */
    cs_real_t  r = 0.5 * _detmat33(n);
    cs_real_t  pi = 4*atan(1.0);

    if (r <= -1)
      theta = _overdim*pi;
    else if (r >= 1)
      theta = 0.;
    else
      theta = _overdim*acos(r);

    // eigenvalues computed should satisfy e1 < e2 < e3
    e3 = tr + 2*p*cos(theta);
    e1 = tr + 2*p*cos(theta + 2*pi*_overdim);
    e2 = 3*tr - e1 -e3; // since tr(m) = e1 + e2 + e3
  }

  /* Debug */
  printf(" --msg-- Computed eigenvalues %5.3e < %5.3e < %5.3e\n", e1, e2, e3);

  /* Return values */
  assert(fabs(e1) > 0);
  *eig_ratio = e3/e1;
  *eig_max = e3;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the convex_hull generated by 3 points.
 *         This corresponds to the computation of the surface of a triangle
 *
 * \param[in]  xv
 * \param[in]  xe
 * \param[in]  xf
 *
 * \return the surface of a triangle
 */
/*----------------------------------------------------------------------------*/

double
cs_surftri(const cs_real_3_t  xv,
           const cs_real_3_t  xe,
           const cs_real_3_t  xf)
{
  int  k;
  cs_real_3_t  u, v, cp;

  double  area = 0.0;

  for (k = 0; k < 3; k++) {
    u[k] = xe[k] - xv[k];
    v[k] = xf[k] - xv[k];
  }
  _cp3(u, v, &cp);
  area = 0.5*_n3(cp);

  return  area;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume of the convex_hull generated by 4 points.
 *         This is equivalent to the computation of the volume of a tetrahedron
 *
 * \param[in]  xv
 * \param[in]  xe
 * \param[in]  xf
 * \param[in]  xc
 *
 * \return the volume of the tetrahedron.
 */
/*----------------------------------------------------------------------------*/

double
cs_voltet(const cs_real_3_t   xv,
          const cs_real_3_t   xe,
          const cs_real_3_t   xf,
          const cs_real_3_t   xc)
{
  double  lev, lef, lec;
  cs_real_3_t  uev, uef, uec, ucp;

  double  vol = 0.0;

  _lenunit3(xe, xv, &lev, &uev);
  _lenunit3(xe, xf, &lef, &uef);
  _lenunit3(xe, xc, &lec, &uec);
  _cp3(uev, uef, &ucp);
  vol = _oversix * lev * lef * lec * fabs(_dp3(ucp, uec));

  return  vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute alpha*x + beta*y = p_z
 *
 * \param[in]     size    vector dimension
 * \param[in]     alpha   coefficient for x vector
 * \param[in]     x       first vector
 * \param[in]     beta    coefficient for y vector
 * \param[in]     y       second vector
 * \param[in,out] p_z     resulting vector (allocated if NULL)
 * \param[in]     reset   reset z vector before computation
 */
/*----------------------------------------------------------------------------*/

void
cs_daxpy(int                size,
         double             alpha,
         const cs_real_t    x[],
         cs_real_t          beta,
         const cs_real_t    y[],
         cs_real_t         *p_z[],
         _Bool              reset)
{
  int  i;

  cs_real_t  *z = *p_z;

  if (size < 1)
    return;

  /* Sanity check */
  assert(x != NULL && y != NULL);

  if (z == NULL)
    reset = true, BFT_MALLOC(z, size, cs_real_t);

  if (reset)
    for (i = 0; i < size; i++) z[i] = 0;

  if (fabs(alpha) < DBL_MIN && fabs(beta) < DBL_MIN)
    return;

  if (fabs(alpha) > DBL_MIN && fabs(beta) > DBL_MIN) {
    for (i = 0; i < size; i++)
      z[i] += alpha*x[i] + beta*y[i];
  }
  else {
    if (fabs(beta) > DBL_MIN)
      for (i = 0; i < size; i++) z[i] += beta*y[i];
    else
      for (i = 0; i < size; i++) z[i] += alpha*y[i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two vectors of dimension "size"
 *         This algorithm tries to reduce round-off error thanks to
 *          intermediate sums.
 *
 *  \param[in] size    vector dimension
 *  \param[in] v       first vector
 *  \param[in] w       second vector
 *
 * \return  the dot product of two vectors
 */
/*----------------------------------------------------------------------------*/

double
cs_dp(int           size,
      const double  v[],
      const double  w[])
{
  int  i, k, sl, test_size, new_size;

  if (size < 1)
    return 0.0;

  if (v == NULL || w == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Array not allocated. Stop dot product computation.\n"));

  /* Initialize sum and index */
  _op_subsum.idx[0] = 0;
  for (i = 0; i < _op_subsum.size; i++) {
    _op_subsum.idx[i+1] = 0;
    _op_subsum.sums[i] = 0.0;
  }

  /* Compute slice size */
  sl = size/_op_subsum.size;
  if (sl % _op_subsum.size != 0)
    sl += 1;
  if (sl == 0)
    sl = 1;

  /* Build index */
  for (i = 0; i < _op_subsum.size; i++) {
    if (_op_subsum.idx[i] < size) {
      _op_subsum.idx[i+1] = _op_subsum.idx[i] + sl;
      if (_op_subsum.idx[i+1] > size) _op_subsum.idx[i+1] = size;
    }
    else
      _op_subsum.idx[i+1] = size;
  }
  _op_subsum.idx[_op_subsum.size] = size;

  /* Compute sums by slice */
  for (k = 0; k < _op_subsum.size; k++)
    for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
      _op_subsum.sums[k] += v[i]*w[i];

  /* Aggregate sub-sums */
  test_size = _op_subsum.size;
  while (test_size > 1) {
    new_size = test_size / 2;
    for (k = 0; k < new_size; k++)
      _op_subsum.sums[k] = _op_subsum.sums[2*k] + _op_subsum.sums[2*k+1];
    if (test_size % 2 != 0)
      _op_subsum.sums[new_size] = _op_subsum.sums[test_size];
    test_size = new_size;
  }

  return _op_subsum.sums[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm 2 of a vector of size len
 *         This algorithm tries to reduce round-off error thanks to
 *          intermediate sums.
 *
 *  \param[in] len     vector dimension
 *  \param[in] v       vector
 *
 * \return  the euclidean norm of a vector
 */
/*----------------------------------------------------------------------------*/

double
cs_euclidean_norm(int            len,
                  const double   v[])
{
  double  n2 = DBL_MAX;

  if (len < 1 || v == NULL)
    return 0.0;

  n2 = cs_dp(len, v, v);
  if (n2 > -DBL_MIN)
    n2 = sqrt(n2);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop norm computation. Norm value is < 0 !\n"));

  return n2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by default the sum of the elements of an array of double
 *         Additional operation are also possible: square, abs
 *         This algorithm tries to reduce round-off errors thanks to
 *         intermediate sums.
 *
 *  \param[in] size    array dimension
 *  \param[in] v       values
 *  \param[in] w       weights (possibly NULL)
 *  \param[in] op      operation to do when doing the sum
 *
 * \return the sum (with possibly additional op) of a vector
 */
/*----------------------------------------------------------------------------*/

double
cs_sum(cs_lnum_t              size,
       const double           v[],
       const double           w[],
       cs_toolbox_type_sum_t  op)
{
  int  i, k, sl, test_size, new_size;

  if (size == 0)
    return 0.0;

  /* Initialize sum and index */
  _op_subsum.idx[0] = 0;
  for (i = 0; i < _op_subsum.size; i++)
    _op_subsum.idx[i+1] = 0, _op_subsum.sums[i] = 0.0;

  /* Compute slice size */
  sl = size/_op_subsum.size;
  if (sl % _op_subsum.size != 0) sl += 1;
  if (sl == 0) sl = 1;

  /* Build index */
  for (i = 0; i < _op_subsum.size; i++) {
    if (_op_subsum.idx[i] < size) {
      _op_subsum.idx[i+1] = _op_subsum.idx[i] + sl;
      if (_op_subsum.idx[i+1] > size)
        _op_subsum.idx[i+1] = size;
    }
    else
      _op_subsum.idx[i+1] = size;
  }
  _op_subsum.idx[_op_subsum.size] = size;

  if (op == CS_TOOLBOX_WSUM    ||
      op == CS_TOOLBOX_WSUMABS ||
      op == CS_TOOLBOX_WSUM2)
    if (w == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" Weighted operation requested but weigths not allocated.\n"
                  " Stop execution.\n"));

  /* Compute sums by slice */
  switch(op) {

  case CS_TOOLBOX_SUM:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += v[i];
    break;

  case CS_TOOLBOX_WSUM:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += w[i]*v[i];
    break;

  case CS_TOOLBOX_SUMABS:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += fabs(v[i]);
    break;

  case CS_TOOLBOX_WSUMABS:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += w[i]*fabs(v[i]);
    break;

  case CS_TOOLBOX_SUM2:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += v[i]*v[i];
    break;

  case CS_TOOLBOX_WSUM2:
    for (k = 0; k < _op_subsum.size; k++)
      for (i = _op_subsum.idx[k]; i < _op_subsum.idx[k+1]; i++)
        _op_subsum.sums[k] += w[i]*v[i]*v[i];
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("  Undefined operation. Stop sum computation.\n"));

  }

  /* Aggregate sub-sums */
  test_size = _op_subsum.size;
  while (test_size > 1) {
    new_size = test_size / 2;
    for (k = 0; k < new_size; k++)
      _op_subsum.sums[k] = _op_subsum.sums[2*k] + _op_subsum.sums[2*k+1];
    if (test_size % 2 != 0)
      _op_subsum.sums[new_size] = _op_subsum.sums[test_size];
    test_size = new_size;
  }

  return _op_subsum.sums[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate or reallocate a temporary buffer structure
 *
 * \param[in]       bufsize   reference size
 * \param[in, out]  p_tb      pointer to the temporary structure to allocate
 */
/*----------------------------------------------------------------------------*/

void
cs_tmpbuf_alloc(size_t           bufsize,
                cs_tmpbuf_t    **p_tb)
{
  cs_tmpbuf_t  *tb = *p_tb;

  if (bufsize == 0)
    return;

  if (tb == NULL) { /* Creation */
    BFT_MALLOC(tb, 1, cs_tmpbuf_t);
    tb->bufsize = bufsize;
    BFT_MALLOC(tb->buf, bufsize, char);
  }
  else { /* Reallocation if needed */
    if (tb->bufsize < bufsize) {
      tb->bufsize = bufsize;
      BFT_REALLOC(tb->buf, bufsize, char);
    }
  }

  /* Returns buffer structure */
  *p_tb = tb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a temporary buffer structure
 *
 * \param[in]  tb      pointer to the temporary structure to free
 *
 * \returns NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_tmpbuf_t *
cs_tmpbuf_free(cs_tmpbuf_t  *tb)
{
  if (tb == NULL)
    return tb;

  BFT_FREE(tb->buf);
  BFT_FREE(tb);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a private structure for this file used for
 *         reducing round-off errors during summation
 *
 *  \param[in] ref_size    reference array dimension
 */
/*----------------------------------------------------------------------------*/

void
cs_toolbox_init(cs_lnum_t      ref_size)
{
  double  invln2 = 1/log(2);
  double  estimate = log(ref_size)*invln2;
  int  power = floor(log(estimate)*invln2);
  int  size = (1 << power);

  /* Compute a number of sub sums according to a reference size */
  _op_subsum.size = CS_MAX(2, size);

  BFT_MALLOC(_op_subsum.idx, _op_subsum.size + 1, int);
  BFT_MALLOC(_op_subsum.sums, _op_subsum.size, double);

  printf("# N_SUB_SUMS      %d\n", _op_subsum.size);
  bft_printf(" -sla- n_subsums      %d\n", _op_subsum.size);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a private structure for this file used for reducing round-off
 *         errors during summation
 */
/*----------------------------------------------------------------------------*/

void
cs_toolbox_finalize(void)
{
  assert(_op_subsum.size > 0);

  BFT_FREE(_op_subsum.idx);
  BFT_FREE(_op_subsum.sums);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_toolbox_locmat_t structure
 *
 * \param[in]  n_max_ent    max number of entities
 *
 * \return  a new allocated cs_toolbox_locmat_t structure
 */
/*----------------------------------------------------------------------------*/

cs_toolbox_locmat_t *
cs_toolbox_locmat_create(int   n_max_ent)
{
  int  i;

  cs_toolbox_locmat_t  *lm = NULL;

  BFT_MALLOC(lm, 1, cs_toolbox_locmat_t);

  lm->n_max_ent = n_max_ent;
  lm->n_ent = 0;
  lm->ids = NULL;
  lm->mat = NULL;

  if (n_max_ent > 0) {

    int  msize = n_max_ent*n_max_ent;

    BFT_MALLOC(lm->ids, n_max_ent, cs_lnum_t);
    for (i = 0; i < n_max_ent; i++)
      lm->ids[i] = 0;

    BFT_MALLOC(lm->mat, msize, double);
    for (i = 0; i < msize; i++)
      lm->mat[i] = 0;

  }

  return lm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local discrete Hodge operator
 *
 * \param[in]    parent_id  id of the related parent entity
 * \param[in]    lm         pointer to the cs_sla_locmat_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_toolbox_locmat_dump(int                         parent_id,
                       const cs_toolbox_locmat_t  *lm)
{
  int  i, j;

  bft_printf("\n  << parent id: %d >>\n", parent_id);

  /* List sub-entity ids */
  for (i = 0; i < lm->n_ent; i++)
    bft_printf(" %9d", lm->ids[i]);
  bft_printf("\n");

  for (i = 0; i < lm->n_ent; i++) {
    bft_printf(" %5d", lm->ids[i]);
    for (j = 0; j < lm->n_ent; j++)
      bft_printf(" % .4e", lm->mat[i*lm->n_ent+j]);
    bft_printf("\n");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_toolbox_locmat_t structure
 *
 * \param[in]  lm    pointer to a cs_toolbox_locmat_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_toolbox_locmat_t *
cs_toolbox_locmat_free(cs_toolbox_locmat_t  *lm)
{
  if (lm == NULL)
    return lm;

  if (lm->n_max_ent > 0) {
    BFT_FREE(lm->ids);
    BFT_FREE(lm->mat);
  }

  BFT_FREE(lm);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute some simple statistics from an array
 *
 * \param[in]  n_elts      number of couples in data
 * \param[in]  stride      size of a couple of data
 * \param[in]  datatype    datatype
 * \param[in]  indata      buffer containing input data
 * \param[in]  do_abs      analyse the absolute value of indata
 *
 * \return a cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

cs_data_info_t
cs_analysis_data(cs_lnum_t       n_elts,
                 int             stride,
                 cs_datatype_t   datatype,
                 const void     *indata,
                 _Bool           do_abs)
{
  int  i, j;

  _Bool  deallocate = false;
  cs_data_info_t  info = _init_dinfo(datatype);

  if (indata == NULL)
    return info;

  switch (datatype) { // Only for double
  case CS_DOUBLE:
    {
      const cs_real_t  *data = (const cs_real_t *)indata;
      cs_real_t  *values = NULL;

      if (stride == 3) { // compute norm

        cs_real_3_t  v;

        BFT_MALLOC(values, n_elts, cs_real_t);
        deallocate = true;
        for (i = 0; i < n_elts; i++) {
          for (j = 0; j < 3; j++)
            v[j] = data[stride*i+j];
          values[i] = _n3(v);
        }
        _compute_info_double(n_elts, values, &info);

      }
      else {  /* stride = 1 */

        assert(stride==1);
        if (do_abs) {
          BFT_MALLOC(values, n_elts, cs_real_t);
          deallocate = true;
          for (i = 0; i < n_elts; i++)
            values[i] = fabs(data[i]);
          _compute_info_double(n_elts, values, &info);
        }
        else
          _compute_info_double(n_elts, data, &info);

      } // stride

      if (deallocate)  BFT_FREE(values);
    }
    break;

  case CS_INT32:
    {
      const cs_lnum_t  *data = (const cs_lnum_t *)indata;
      cs_lnum_t  *numbers = NULL;

      if (do_abs) {
        BFT_MALLOC(numbers, n_elts, cs_lnum_t);
        for (i = 0; i < n_elts; i++)
          numbers[i] = CS_ABS(data[i]);
        _compute_info_int32(n_elts, numbers, &info);
        BFT_FREE(numbers);
      }
      else
        _compute_info_int32(n_elts, data, &info);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_data_info_t structure
 *
 * \param[in]  name        filename if not NULL
 * \param[in]  f           output file if not NULL
 * \param[in]  n_elts      number of couples in data
 * \param[in]  datatype    datatype
 * \param[in]  dinfo       cs_data_info_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_data_info_dump(const char              *name,
                  FILE                    *f,
                  cs_lnum_t                n_elts,
                  cs_datatype_t            datatype,
                  const cs_data_info_t     dinfo)
{
  FILE  *_f = f;
  _Bool close_file = false;

  if (f == NULL) {
    if (name == NULL) _f = stdout;
    else
      _f = fopen(name, "w"), close_file = true;
  }

  fprintf(_f, "\n");
  if (name != NULL)
    fprintf(_f, " -dat- name          %s\n", name);
  fprintf(_f, " -dat- n_elts        %d\n", n_elts);

  switch (datatype) {
  case CS_DOUBLE:
    fprintf(_f, " -dat- minimum    %- 9.6e\n", dinfo.min.value);
    fprintf(_f, " -dat- maximum    %- 9.6e\n", dinfo.max.value);
    break;
  case CS_INT32:
    fprintf(_f, " -dat- minimum    %10d\n", dinfo.min.number);
    fprintf(_f, " -dat- maximum    %10d\n", dinfo.max.number);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid datatype for analysing data.\n"));
  }

  fprintf(_f, " -dat- mean            %- 9.6e\n", dinfo.mean);
  fprintf(_f, " -dat- sigma           %- 9.6e\n", dinfo.sigma);
  fprintf(_f, " -dat- euclidean norm  %- 9.6e\n", dinfo.euclidean_norm);
  fprintf(_f, "\n");

  fflush(_f);
  if (close_file)  fclose(_f);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
