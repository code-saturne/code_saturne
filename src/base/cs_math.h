#ifndef __CS_MATH_H__
#define __CS_MATH_H__

/*============================================================================
 * Mathematical base functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if  (defined(__NVCC__) && defined(__CUDA_ARCH__)) \
  || defined(SYCL_LANGUAGE_VERSION) \
  || defined(HAVE_OPENMP_TARGET)
#include <float.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
#include "bft/bft_error.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/* Symmetric tensor component name */

typedef enum {

  XX,
  YY,
  ZZ,
  XY,
  YZ,
  XZ

} cs_math_sym_tensor_component_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Numerical constants */

#if  (defined(__NVCC__) && defined(__CUDA_ARCH__)) \
  || defined(SYCL_LANGUAGE_VERSION) \
  || defined(HAVE_OPENMP_TARGET)

/* On GPU, global variables are usually not accessible. */

#define cs_math_zero_threshold FLT_MIN
#define cs_math_epzero 1e-12
#define cs_math_infinite_r 1.e30
#define cs_math_big_r 1.e12
#define cs_math_pi 3.14159265358979323846

#else

/* General constants accessible on CPU */

extern const cs_real_t cs_math_zero_threshold;
extern const cs_real_t cs_math_epzero;
extern const cs_real_t cs_math_infinite_r;
extern const cs_real_t cs_math_big_r;
extern const cs_real_t cs_math_pi;

#endif

#if !(defined(__NVCC__) && defined(__CUDA_ARCH__))

extern const cs_real_t cs_math_1ov3;
extern const cs_real_t cs_math_2ov3;
extern const cs_real_t cs_math_4ov3;
extern const cs_real_t cs_math_5ov3;
extern const cs_real_t cs_math_1ov6;
extern const cs_real_t cs_math_1ov12;
extern const cs_real_t cs_math_1ov24;

/* Identity matrix in dimension 3 */
static const cs_real_33_t cs_math_33_identity = {{1., 0., 0.,},
                                                 {0., 1., 0.},
                                                 {0., 0., 1.}};
static const cs_real_6_t cs_math_sym_33_identity  = {1., 1., 1., 0. ,0., 0.};

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*=============================================================================
 * Templated inline functions
 *============================================================================*/

#if defined(__cplusplus)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute \f$ (\vect{x}_b - \vect{x}_a) \cdot \vect{x}_c \f$
 *
 * \tparam T input value type
 * \tparam U input value type
 *
 * \param[in]  xa   first coordinate
 * \param[in]  xb   second coordinate
 * \param[in]  xc   third coordinate
 *
 * \return the dot product
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_distance_dot_product(const T  xa[3],
                               const T  xb[3],
                               const U  xc[3])
{
  return ((xb[0] - xa[0])*xc[0]+(xb[1] - xa[1])*xc[1]+(xb[2] - xa[2])*xc[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of two vectors of 3 real values.
 *
 * \tparam T input value type
 * \tparam U input value type
 *
 * \param[in]   u   vector of 3 real values
 * \param[in]   v   vector of 3 real values
 *
 * \return the resulting dot product u.v.
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_dot_product(const T  u[3],
                      const U  v[3])
{
  double uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

  return uv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm of a vector of dimension 3
 *
 * \tparam T  input value type
 *
 * \param[in]  v
 *
 * \return the value of the norm
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_norm(const T  v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of a symmetric tensor t with two vectors,
 *        n1 and n2.
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]     n1    vector of 3 real values
 * \param[in]     t     tensor of 6 real values
 *                      [ 0 3 5 ]
 *                      [ 3 1 4 ]
 *                      [ 5 4 2 ]
 * \param[in]     n2    vector of 3 real values
 *
 * \return the resulting dot product n1.t.n2.
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_sym_33_3_dot_product(const T  n1[3],
                               const U  t[6],
                               const V  n2[3])
{
  return (  n1[0] * (t[0]*n2[0] + t[3]*n2[1] + t[5]*n2[2])
          + n1[1] * (t[3]*n2[0] + t[1]*n2[1] + t[4]*n2[2])
          + n1[2] * (t[5]*n2[0] + t[4]*n2[1] + t[2]*n2[2]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the square norm of a vector of 3 real values.
 *
 * \tparam T  input value type
 *
 * \param[in]     v             vector of 3 real values
 *
 * \return square norm of v.
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_square_norm(const T v[3])
{
  cs_real_t v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

  return v2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalise a vector of 3 real values.
 *
 * \tparam T input value type
 * \tparam U input value type
 *
 * To normalize in-place, \p vin and \p vout may point to the same array.
 *
 * \param[in]     vin           vector
 * \param[out]    vout          normalized vector
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U>
CS_F_HOST_DEVICE static inline void
cs_math_3_normalize(const T  vin[3],
                    U        vout[3])
{
  cs_real_t norm = cs_math_3_norm(vin);

  cs_real_t inv_norm = ((norm > cs_math_zero_threshold) ?  1. / norm : 0);

  vout[0] = inv_norm * vin[0];
  vout[1] = inv_norm * vin[1];
  vout[2] = inv_norm * vin[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a symmetric matrix of 3x3 real values by
 * a vector of 3 real values.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE inline void
cs_math_sym_33_3_product(const T       m[6],
                         const U       v[3],
                         V  *restrict  mv)
{
  mv[0] = m[0]*v[0] + m[3]*v[1] + m[5]*v[2];
  mv[1] = m[3]*v[0] + m[1]*v[1] + m[4]*v[2];
  mv[2] = m[5]*v[0] + m[4]*v[1] + m[2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the dot product with a normal vector to the normal direction
 *        to a vector.
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]       n       normalised face normal vector
 * \param[in]       factor  factor
 * \param[in, out]  v       vector to be scaled
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE inline void
cs_math_3_normal_scaling(const T  n[3],
                         U        factor,
                         V        v[3])
{
  cs_real_t v_dot_n = (factor -1.) * cs_math_3_dot_product(v, n);
  for (int i = 0; i < 3; i++)
    v[i] += v_dot_n * n[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of a tensor t with two vectors, n1 and n2.
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]     n1    vector of 3 real values
 * \param[in]     t     tensor of 3x3 real values
 * \param[in]     n2    vector of 3 real values
 *
 * \return the resulting dot product n1.t.n2.
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE inline cs_real_t
cs_math_3_33_3_dot_product(const T  n1[3],
                           const U  t[3][3],
                           const V  n2[3])
{
  cs_real_t n_t_n
    = (  n1[0]*t[0][0]*n2[0] + n1[1]*t[1][0]*n2[0] + n1[2]*t[2][0]*n2[0]
       + n1[0]*t[0][1]*n2[1] + n1[1]*t[1][1]*n2[1] + n1[2]*t[2][1]*n2[1]
       + n1[0]*t[0][2]*n2[2] + n1[1]*t[1][2]*n2[2] + n1[2]*t[2][2]*n2[2]);
  return n_t_n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Orthogonal projection of a vector with respect to a normalised
 *        vector.
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]   n     normal vector direction
 * \param[in]   v     vector to be projected
 * \param[out]  vout  projection
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE inline void
cs_math_3_orthogonal_projection(const T      n[3],
                                const U      v[3],
                                V  *restrict vout)
{
  vout[0] =  v[0]*(1.-n[0]*n[0])- v[1]*    n[1]*n[0] - v[2]*    n[2]*n[0];
  vout[1] = -v[0]*    n[0]*n[1] + v[1]*(1.-n[1]*n[1])- v[2]*    n[2]*n[1];
  vout[2] = -v[0]*    n[0]*n[2] - v[1]*    n[1]*n[2] + v[2]*(1.-n[2]*n[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the cross product of two vectors of 3 real values.
 *
 * \tparam T input value type
 * \tparam U input value type
 * \tparam V input value type
 *
 * \param[in]     u    vector of 3 real values
 * \param[in]     v    vector of 3 real values
 * \param[out]    uv   cross-product of u an v
 */
/*----------------------------------------------------------------------------*/

template <typename T, typename U, typename V>
CS_F_HOST_DEVICE static inline void
cs_math_3_cross_product(const T       u[3],
                        const U       v[3],
                        V  *restrict  uv)
{
  uv[0] = u[1]*v[2] - u[2]*v[1];
  uv[1] = u[2]*v[0] - u[0]*v[2];
  uv[2] = u[0]*v[1] - u[1]*v[0];
}

namespace cs {

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return absolute value of a number.
 *
 * Specialized versions for floating point numbers should be faster than
 * base version from cs_defs.h.
 *
 * \param[in]  a  value
 *
 * \return absolute value of a
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline double
abs(const double  a)
{
  return fabs(a);
}

CS_F_HOST_DEVICE inline float
abs(const float  a)
{
  return fabsf(a);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return minimum of two values.
 *
 * Specialized versions for floating point numbers should be faster than
 * base version from cs_defs.h.
 *
 * \param[in]  a   first value
 * \param[in]  b   second value
 *
 * \return minimum of a and b
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline double
min(const double  a,
    const double  b)
{
  return fmin(a, b);
}

CS_F_HOST_DEVICE inline float
min(const float  a,
    const float  b)
{
  return fminf(a, b);
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return maximum of two values.
 *
 * Specialized versions for floating point numbers should be faster than
 * base version from cs_defs.h.
 *
 * \param[in]  a   first value
 * \param[in]  b   second value
 *
 * \return maximum of a and b
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline double
max(const double  a,
    const double  b)
{
  return fmax(a, b);
}

CS_F_HOST_DEVICE inline float
max(const float  a,
    const float  b)
{
  return fmaxf(a, b);
}

} // namespace cs

#endif // defined(__cplusplus)

BEGIN_C_DECLS

/*=============================================================================
 * Inline static functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Computes the binomial coefficient of n and k
 *
 * \param[in] n      first argument
 * \param[in] k      second argument
 *
 * \return the value of binom(n)(k)
 */
/*----------------------------------------------------------------------------*/

static inline int
cs_math_binom(int  n,
              int  k)
{
  int ret = 1;
  assert(n >= k);

  const int n_iter = (k > n-k) ? n-k : k;
  for (int j = 1; j <= n_iter; j++, n--) {
    if (n % j == 0)
      ret *= n/j;
    else if (ret % j == 0)
      ret = ret/j*n;
    else
      ret = (ret*n)/j;
  }

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the absolute value of a real value.
 *
 * \param[in]  x  value
 *
 * \return absolute value of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_fabs(cs_real_t  x)
{
  cs_real_t ret = (x <  0) ? -x : x;

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the min value of two real values.
 *
 * \param[in]  x, y  values
 *
 * \return the min value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_fmin(cs_real_t  x,
             cs_real_t  y)
{
  cs_real_t ret = (x <  y) ? x : y;

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the max value of two real values.
 *
 * \param[in]  x, y  values
 *
 * \return the max value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_fmax(cs_real_t  x,
             cs_real_t  y)
{
  cs_real_t ret = (x <  y) ? y : x;

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clamp function for a given scalar value.
 *
 * \param[in] x    initial value
 * \param[in] xmin min value for clamping
 * \param[in] xmax max value for clamping
 *
 * \return clamped value which is 'x' if xmin < x < xmax or lower/upper bound
 * otherwise
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_clamp(cs_real_t x,
              cs_real_t xmin,
              cs_real_t xmax)
{
  cs_real_t ret = cs_math_fmin(xmax, cs_math_fmax(xmin, x));

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square of a real value.
 *
 * \param[in]  x  value
 *
 * \return the square of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_sq(cs_real_t  x)
{
  return x*x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square of a real value.
 *
 * \param[in]  x  value
 *
 * \return the square of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_pow2(cs_real_t  x)
{
  return x*x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cube of a real value.
 *
 * \param[in]  x  value
 *
 * \return the cube of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_pow3(cs_real_t  x)
{
  return x*x*x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the 4-th power of a real value.
 *
 * \param[in]  x  value
 *
 * \return the 4th power of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_pow4(cs_real_t  x)
{
  cs_real_t x2 = x * x;
  return x2 * x2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the 5-th power of a real value.
 *
 * \param[in]  x  value
 *
 * \return the 5th power of the given value
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_pow5(cs_real_t  x)
{
  cs_real_t x2 = x * x;
  return x * x2 * x2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the (euclidean) distance between two points xa and xb in
 *         a Cartesian coordinate system of dimension 3.
 *
 * \param[in]  xa   first coordinate
 * \param[in]  xb   second coordinate
 *
 * \return the length between two points xa and xb
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_distance(const cs_real_t  xa[3],
                   const cs_real_t  xb[3])
{
  cs_real_t  v[3];

  v[0] = xb[0] - xa[0];
  v[1] = xb[1] - xa[1];
  v[2] = xb[2] - xa[2];

  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute \f$ (\vect{x}_b - \vect{x}_a) \cdot \vect{x}_c \f$
 *
 * \param[in]  xa   first coordinate
 * \param[in]  xb   second coordinate
 * \param[in]  xc   third coordinate
 *
 * \return the dot product
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_distance_dot_product(const cs_real_t  xa[3],
                               const cs_real_t  xb[3],
                               const cs_real_t  xc[3])
{
  return ((xb[0] - xa[0])*xc[0]+(xb[1] - xa[1])*xc[1]+(xb[2] - xa[2])*xc[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the squared distance between two points xa and xb in
 *         a Cartesian coordinate system of dimension 3.
 *
 * \param[in]  xa   first coordinate
 * \param[in]  xb   second coordinate
 *
 * \return the square length between two points xa and xb
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_square_distance(const cs_real_t  xa[3],
                          const cs_real_t  xb[3])
{
  cs_real_t v[3] = {xb[0] - xa[0],
                    xb[1] - xa[1],
                    xb[2] - xa[2]};

  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of two vectors of 3 real values.
 *
 * \param[in]   u   vector of 3 real values
 * \param[in]   v   vector of 3 real values
 *
 * \return the resulting dot product u.v.
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_dot_product(const cs_real_t  u[3],
                      const cs_real_t  v[3])
{
  cs_real_t uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

  return uv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of a tensor t with two vectors, n1 and n2.
 *
 * \param[in]     n1    vector of 3 real values
 * \param[in]     t     tensor of 3x3 real values
 * \param[in]     n2    vector of 3 real values
 *
 * \return the resulting dot product n1.t.n2.
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_33_3_dot_product(const cs_real_t  n1[3],
                           const cs_real_t  t[3][3],
                           const cs_real_t  n2[3])
{
  cs_real_t n_t_n
    = (  n1[0]*t[0][0]*n2[0] + n1[1]*t[1][0]*n2[0] + n1[2]*t[2][0]*n2[0]
       + n1[0]*t[0][1]*n2[1] + n1[1]*t[1][1]*n2[1] + n1[2]*t[2][1]*n2[1]
       + n1[0]*t[0][2]*n2[2] + n1[1]*t[1][2]*n2[2] + n1[2]*t[2][2]*n2[2]);
  return n_t_n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of a symmetric tensor t with two vectors,
 *        n1 and n2.
 *
 * \param[in]     n1    vector of 3 real values
 * \param[in]     t     tensor of 6 real values
 *                      [ 0 3 5 ]
 *                      [ 3 1 4 ]
 *                      [ 5 4 2 ]
 * \param[in]     n2    vector of 3 real values
 *
 * \return the resulting dot product n1.t.n2.
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_sym_33_3_dot_product(const cs_real_t  n1[3],
                               const cs_real_t  t[6],
                               const cs_real_t  n2[3])
{
  return (  n1[0] * (t[0]*n2[0] + t[3]*n2[1] + t[5]*n2[2])
          + n1[1] * (t[3]*n2[0] + t[1]*n2[1] + t[4]*n2[2])
          + n1[2] * (t[5]*n2[0] + t[4]*n2[1] + t[2]*n2[2]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm of a vector of dimension 3
 *
 * \param[in]  v
 *
 * \return the value of the norm
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_norm(const cs_real_t  v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the square norm of a vector of 3 real values.
 *
 * \param[in]     v             vector of 3 real values
 *
 * \return square norm of v.
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_square_norm(const cs_real_t v[3])
{
  cs_real_t v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

  return v2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalise a vector of 3 real values.
 *
 * To normalize in-place, \p vin and \p vout may point to the same array.
 *
 * \param[in]     vin           vector
 * \param[out]    vout          normalized vector
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_normalize(const cs_real_t  vin[3],
                    cs_real_t        vout[3])
{
  cs_real_t norm = cs_math_3_norm(vin);

  cs_real_t inv_norm = ((norm > cs_math_zero_threshold) ?  1. / norm : 0);

  vout[0] = inv_norm * vin[0];
  vout[1] = inv_norm * vin[1];
  vout[2] = inv_norm * vin[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalise a vector of 3 real values and clip the norm using a
 * threshold value.
 *
 * To normalize in-place, \p vin and \p vout may point to the same array.
 *
 * \param[in]     vin           vector
 * \param[in]     thres         threshold for normalization
 * \param[out]    vout          normalized vector
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_normalize_threshold(const cs_real_t vin[3],
                              const cs_real_t thres,
                              cs_real_t       vout[3])
{
  cs_real_t norm = cs_math_3_norm(vin);

  cs_real_t inv_norm = ((norm > thres) ? 1. / norm : 1. / thres);

  vout[0] = inv_norm * vin[0];
  vout[1] = inv_norm * vin[1];
  vout[2] = inv_norm * vin[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Orthogonal projection of a vector with respect to a normalised
 *        vector.
 *
 * \param[in]   n     normal vector direction
 * \param[in]   v     vector to be projected
 * \param[out]  vout  projection
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_orthogonal_projection(const cs_real_t      n[3],
                                const cs_real_t      v[3],
                                cs_real_t  *restrict vout)
{
  vout[0] =  v[0]*(1.-n[0]*n[0])- v[1]*    n[1]*n[0] - v[2]*    n[2]*n[0];
  vout[1] = -v[0]*    n[0]*n[1] + v[1]*(1.-n[1]*n[1])- v[2]*    n[2]*n[1];
  vout[2] = -v[0]*    n[0]*n[2] - v[1]*    n[1]*n[2] + v[2]*(1.-n[2]*n[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the dot product with a normal vector to the normal direction
 *        to a vector.
 *
 * \param[in]       n       normalised face normal vector
 * \param[in]       factor  factor
 * \param[in, out]  v       vector to be scaled
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_normal_scaling(const cs_real_t  n[3],
                         cs_real_t        factor,
                         cs_real_t        v[3])
{
  cs_real_t v_dot_n = (factor -1.) * cs_math_3_dot_product(v, n);
  for (int i = 0; i < 3; i++)
    v[i] += v_dot_n * n[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the dot product with a normal vector to the normal,normal
 * component of a tensor:
 * t += factor * n.t.n n(x)n
 *
 * \param[in]     n             normalised face normal vector
 * \param[in]     factor        factor
 * \param[in,out] t             matrix to be scaled
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_normal_scaling_add(const cs_real_t  n[3],
                              cs_real_t        factor,
                              cs_real_t        t[3][3])
{
  cs_real_t n_t_n = (factor -1.) *
    ( n[0] * t[0][0] * n[0] + n[1] * t[1][0] * n[0] + n[2] * t[2][0] * n[0]
    + n[0] * t[0][1] * n[1] + n[1] * t[1][1] * n[1] + n[2] * t[2][1] * n[1]
    + n[0] * t[0][2] * n[2] + n[1] * t[1][2] * n[2] + n[2] * t[2][2] * n[2]);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      t[i][j] += n_t_n * n[i] * n[j];
  }
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a matrix of 3x3 real values by a vector of 3
 * real values.
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_3_product(const cs_real_t      m[3][3],
                     const cs_real_t      v[3],
                     cs_real_t  *restrict mv)
{
  mv[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  mv[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  mv[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a matrix of 3x3 real values by a vector of 3
 * real values add.
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[in,out] mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_3_product_add(const cs_real_t       m[3][3],
                         const cs_real_t       v[3],
                         cs_real_t  *restrict  mv)
{
  mv[0] += m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  mv[1] += m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  mv[2] += m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of the transpose of a matrix of 3x3 real
 * values by a vector of 3 real values.
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33t_3_product(const cs_real_t       m[3][3],
                      const cs_real_t       v[3],
                      cs_real_t  *restrict  mv)
{
  mv[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  mv[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  mv[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a symmetric matrix of 3x3 real values by
 * a vector of 3 real values.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_3_product(const cs_real_t       m[6],
                         const cs_real_t       v[3],
                         cs_real_t  *restrict  mv)
{
  mv[0] = m[0]*v[0] + m[3]*v[1] + m[5]*v[2];
  mv[1] = m[3]*v[0] + m[1]*v[1] + m[4]*v[2];
  mv[2] = m[5]*v[0] + m[4]*v[1] + m[2]*v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a symmetric matrix of 3x3 real values by
 * a vector of 3 real values and add it to the vector.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_3_product_add(const cs_real_t       m[6],
                             const cs_real_t       v[3],
                             cs_real_t  *restrict  mv)
{
  mv[0] += m[0] * v[0] + m[3] * v[1] + m[5] * v[2];
  mv[1] += m[3] * v[0] + m[1] * v[1] + m[4] * v[2];
  mv[2] += m[5] * v[0] + m[4] * v[1] + m[2] * v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of two symmetric matrices of 3x3 real values
 * and take the trace.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 * \return trace
 *
 * \param[in]     m1            matrix of 3x3 real values
 * \param[in]     m2            matrix of 3x3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_sym_33_sym_33_product_trace(const cs_real_t  m1[6],
                                    const cs_real_t  m2[6])
{
  return m1[0]*m2[0] + 2.*m1[3]*m2[3] + 2.*m1[5]*m2[5]
                     +    m1[1]*m2[1] + 2.*m1[4]*m2[4]
                                      +    m1[2]*m2[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the trace of a 3x3 tensor.
 *
 * \param[in]   t   tensor of 3x3 real values
 *
 * \return trace (t[0][0] + t[1][1] + t[2][2])
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_33_trace(const cs_real_t  t[3][3])
{
  return (t[0][0] + t[1][1] + t[2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the trace of a symmetric tensor.
 *
 * \param[in]   t   vector of 6 real values (symmetric tensor)
 *
 * \return trace (t[0] + t[1] + t[2])
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_6_trace(const cs_real_t  t[6])
{
  return (t[0] + t[1] + t[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a matrix of 6x6 real values by
 * a vector of 6 real values.
 *
 * \param[in]     m             matrix of 6x6 real values
 * \param[in]     v             vector of 6 real values
 * \param[out]    mv            vector of 6 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_66_6_product(const cs_real_t       m[6][6],
                     const cs_real_t       v[6],
                     cs_real_t  *restrict  mv)
{
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++)
      mv[i] = m[i][j] * v[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a matrix of 6x6 real values by
 * a vector of 6 real values and add it to the vector.
 *
 * \param[in]     m             matrix of 6x6 real values
 * \param[in]     v             vector of 6 real values
 * \param[out]    mv            vector of 6 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_66_6_product_add(const cs_real_t       m[6][6],
                         const cs_real_t       v[6],
                         cs_real_t  *restrict  mv)
{
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++)
      mv[i] += m[i][j] * v[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the determinant of a 3x3 matrix
 *
 * \param[in]  m    3x3 matrix
 *
 * \return the determinant
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_33_determinant(const cs_real_t   m[3][3])
{
  const cs_real_t  com0 = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  const cs_real_t  com1 = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  const cs_real_t  com2 = m[0][1]*m[1][2] - m[1][1]*m[0][2];

  return m[0][0]*com0 + m[1][0]*com1 + m[2][0]*com2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the determinant of a 3x3 symmetric matrix
 *
 * \param[in]  m    3x3 symmetric matrix
 *
 * \return the determinant
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_sym_33_determinant(const cs_real_t  m[6])
{
  const cs_real_t  com0 = m[1]*m[2] - m[4]*m[4];
  const cs_real_t  com1 = m[4]*m[5] - m[3]*m[2];
  const cs_real_t  com2 = m[3]*m[4] - m[1]*m[5];

  return m[0]*com0 + m[3]*com1 + m[5]*com2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the average of two vector of dimension 3
 *
 * \param[in]     u    vector of 3 real values
 * \param[in]     v    vector of 3 real values
 * \param[out]    uv   average of u an v := (u+v)/2
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_average(const cs_real_t       u[3],
                  const cs_real_t       v[3],
                  cs_real_t  *restrict  uv)
{
  uv[0] = (u[0] + v[0]) / 2.0;
  uv[1] = (u[1] + v[1]) / 2.0;
  uv[2] = (u[2] + v[2]) / 2.0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the cross product of two vectors of 3 real values.
 *
 * \param[in]     u    vector of 3 real values
 * \param[in]     v    vector of 3 real values
 * \param[out]    uv   cross-product of u an v
 */
/*----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 0 /* Bug with O1 or above with icc 15.0.1 20141023 */
#endif

CS_F_HOST_DEVICE static inline void
cs_math_3_cross_product(const cs_real_t       u[3],
                        const cs_real_t       v[3],
                        cs_real_t  *restrict  uv)
{
  uv[0] = u[1]*v[2] - u[2]*v[1];
  uv[1] = u[2]*v[0] - u[0]*v[2];
  uv[2] = u[0]*v[1] - u[1]*v[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the triple product
 *
 * \param[in]     u   vector of 3 real values
 * \param[in]     v   vector of 3 real values
 * \param[in]     w   vector of 3 real values
 *
 * \return the scalar triple product: (uxv).w
 */
/*----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 0 /* Bug with O1 or above with icc 15.0.1 20141023 */
#endif

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_3_triple_product(const cs_real_t  u[3],
                         const cs_real_t  v[3],
                         const cs_real_t  w[3])
{
  return    (u[1]*v[2] - u[2]*v[1]) * w[0]
          + (u[2]*v[0] - u[0]*v[2]) * w[1]
          + (u[0]*v[1] - u[1]*v[0]) * w[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build an orthonormal basis based on a first vector "vect".
 *        axes[0] is vect normalized, while (axes[0], axes[1], axes[23])
 *        is an orthonormal base.
 *
 * \param[in]  vect Vector used to build the orthonormal basis
 * \param[out] axes axes basis (cs_real_t[3][3])
 *
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_3_orthonormal_basis(const cs_real_t vect[3],
                            cs_real_t       axes[3][3])
{
  assert(cs_math_3_norm(vect) > cs_math_zero_threshold);

  // Compute normal - third axis
  cs_math_3_normalize(vect, axes[2]);

  // Compute base's first axis
  // First test projection of Ox
  cs_real_t Ox[3] = {1., 0., 0.};
  cs_real_t w[3] = {0.};

  cs_math_3_orthogonal_projection(axes[2], Ox, w);

  // If Ox projection is null, project Oy
  if (cs_math_3_norm(w) < cs_math_zero_threshold) {
    cs_real_t Oy[3] = {0., 1., 0.};
    cs_math_3_orthogonal_projection(axes[2], Oy, w);
  }

  cs_math_3_normalize(w, axes[0]);

  // Compute base's second axis using cross product
  cs_math_3_cross_product(axes[2], axes[0], axes[1]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 matrix
 *
 * \param[in]  in    matrix to inverse
 * \param[out] out   inversed matrix
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_inv_cramer(const cs_real_t  in[3][3],
                      cs_real_t        out[3][3])
{
  out[0][0] = in[1][1]*in[2][2] - in[2][1]*in[1][2];
  out[0][1] = in[2][1]*in[0][2] - in[0][1]*in[2][2];
  out[0][2] = in[0][1]*in[1][2] - in[1][1]*in[0][2];

  out[1][0] = in[2][0]*in[1][2] - in[1][0]*in[2][2];
  out[1][1] = in[0][0]*in[2][2] - in[2][0]*in[0][2];
  out[1][2] = in[1][0]*in[0][2] - in[0][0]*in[1][2];

  out[2][0] = in[1][0]*in[2][1] - in[2][0]*in[1][1];
  out[2][1] = in[2][0]*in[0][1] - in[0][0]*in[2][1];
  out[2][2] = in[0][0]*in[1][1] - in[1][0]*in[0][1];

  const double  det = in[0][0]*out[0][0]+in[1][0]*out[0][1]+in[2][0]*out[0][2];
  const double  invdet = 1./det;

  out[0][0] *= invdet, out[0][1] *= invdet, out[0][2] *= invdet;
  out[1][0] *= invdet, out[1][1] *= invdet, out[1][2] *= invdet;
  out[2][0] *= invdet, out[2][1] *= invdet, out[2][2] *= invdet;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 matrix in place, using Cramer's rule
 *
 * \param[in, out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_inv_cramer_in_place(cs_real_t  a[3][3])
{
  cs_real_t a00 = a[1][1]*a[2][2] - a[2][1]*a[1][2];
  cs_real_t a01 = a[2][1]*a[0][2] - a[0][1]*a[2][2];
  cs_real_t a02 = a[0][1]*a[1][2] - a[1][1]*a[0][2];
  cs_real_t a10 = a[2][0]*a[1][2] - a[1][0]*a[2][2];
  cs_real_t a11 = a[0][0]*a[2][2] - a[2][0]*a[0][2];
  cs_real_t a12 = a[1][0]*a[0][2] - a[0][0]*a[1][2];
  cs_real_t a20 = a[1][0]*a[2][1] - a[2][0]*a[1][1];
  cs_real_t a21 = a[2][0]*a[0][1] - a[0][0]*a[2][1];
  cs_real_t a22 = a[0][0]*a[1][1] - a[1][0]*a[0][1];

  double det_inv = 1. / (a[0][0]*a00 + a[1][0]*a01 + a[2][0]*a02);

  a[0][0] = a00 * det_inv;
  a[0][1] = a01 * det_inv;
  a[0][2] = a02 * det_inv;
  a[1][0] = a10 * det_inv;
  a[1][1] = a11 * det_inv;
  a[1][2] = a12 * det_inv;
  a[2][0] = a20 * det_inv;
  a[2][1] = a21 * det_inv;
  a[2][2] = a22 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 symmetric matrix (with non-symmetric storage)
 *         in place, using Cramer's rule
 *
 * \param[in, out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_inv_cramer_sym_in_place(cs_real_t  a[3][3])
{
  cs_real_t a00 = a[1][1]*a[2][2] - a[2][1]*a[1][2];
  cs_real_t a01 = a[2][1]*a[0][2] - a[0][1]*a[2][2];
  cs_real_t a02 = a[0][1]*a[1][2] - a[1][1]*a[0][2];
  cs_real_t a11 = a[0][0]*a[2][2] - a[2][0]*a[0][2];
  cs_real_t a12 = a[1][0]*a[0][2] - a[0][0]*a[1][2];
  cs_real_t a22 = a[0][0]*a[1][1] - a[1][0]*a[0][1];

  double det_inv = 1. / (a[0][0]*a00 + a[1][0]*a01 + a[2][0]*a02);

  a[0][0] = a00 * det_inv;
  a[0][1] = a01 * det_inv;
  a[0][2] = a02 * det_inv;
  a[1][0] = a01 * det_inv;
  a[1][1] = a11 * det_inv;
  a[1][2] = a12 * det_inv;
  a[2][0] = a02 * det_inv;
  a[2][1] = a12 * det_inv;
  a[2][2] = a22 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of a symmetric matrix using Cramer's rule.
 *
 * \remark Symmetric matrix coefficients are stored as follows:
 *         (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s      symmetric matrix
 * \param[out]    sout   sout = 1/s1
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_inv_cramer(const cs_real_t       s[6],
                          cs_real_t  *restrict  sout)
{
  double detinv;

  sout[0] = s[1]*s[2] - s[4]*s[4];
  sout[1] = s[0]*s[2] - s[5]*s[5];
  sout[2] = s[0]*s[1] - s[3]*s[3];
  sout[3] = s[4]*s[5] - s[3]*s[2];
  sout[4] = s[3]*s[5] - s[0]*s[4];
  sout[5] = s[3]*s[4] - s[1]*s[5];

  detinv = 1. / (s[0]*sout[0] + s[3]*sout[3] + s[5]*sout[5]);

  sout[0] *= detinv;
  sout[1] *= detinv;
  sout[2] *= detinv;
  sout[3] *= detinv;
  sout[4] *= detinv;
  sout[5] *= detinv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of two 3x3 real valued matrices.
 *
 * \param[in]     m1     matrix of 3x3 real values
 * \param[in]     m2     matrix of 3x3 real values
 * \param[out]    mout   m1.m2 product
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_product(const cs_real_t  m1[3][3],
                   const cs_real_t  m2[3][3],
                   cs_real_t        mout[3][3])
{
  mout[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
  mout[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1];
  mout[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2];

  mout[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
  mout[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1];
  mout[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2];

  mout[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
  mout[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1];
  mout[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute transformation from relative to absolute
 *        reference frame Q^t M Q
 *
 * \param[in]     m      matrix of 3x3 real values
 * \param[in]     q      transformation matrix of 3x3 real values
 * \param[out]    mout   Q^t M Q
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_transform_r_to_a(const cs_real_t  m[3][3],
                            const cs_real_t  q[3][3],
                            cs_real_t        mout[3][3])
{
  /* _m = M.Q */
  cs_real_33_t _m;
  _m[0][0] = m[0][0]*q[0][0] + m[0][1]*q[1][0] + m[0][2]*q[2][0];
  _m[0][1] = m[0][0]*q[0][1] + m[0][1]*q[1][1] + m[0][2]*q[2][1];
  _m[0][2] = m[0][0]*q[0][2] + m[0][1]*q[1][2] + m[0][2]*q[2][2];

  _m[1][0] = m[1][0]*q[0][0] + m[1][1]*q[1][0] + m[1][2]*q[2][0];
  _m[1][1] = m[1][0]*q[0][1] + m[1][1]*q[1][1] + m[1][2]*q[2][1];
  _m[1][2] = m[1][0]*q[0][2] + m[1][1]*q[1][2] + m[1][2]*q[2][2];

  _m[2][0] = m[2][0]*q[0][0] + m[2][1]*q[1][0] + m[2][2]*q[2][0];
  _m[2][1] = m[2][0]*q[0][1] + m[2][1]*q[1][1] + m[2][2]*q[2][1];
  _m[2][2] = m[2][0]*q[0][2] + m[2][1]*q[1][2] + m[2][2]*q[2][2];

  /* mout = Q^t _m */
  mout[0][0] = q[0][0]*_m[0][0] + q[1][0]*_m[1][0] + q[2][0]*_m[2][0];
  mout[0][1] = q[0][0]*_m[0][1] + q[1][0]*_m[1][1] + q[2][0]*_m[2][1];
  mout[0][2] = q[0][0]*_m[0][2] + q[1][0]*_m[1][2] + q[2][0]*_m[2][2];

  mout[1][0] = q[0][1]*_m[0][0] + q[1][1]*_m[1][0] + q[2][1]*_m[2][0];
  mout[1][1] = q[0][1]*_m[0][1] + q[1][1]*_m[1][1] + q[2][1]*_m[2][1];
  mout[1][2] = q[0][1]*_m[0][2] + q[1][1]*_m[1][2] + q[2][1]*_m[2][2];

  mout[2][0] = q[0][2]*_m[0][0] + q[1][2]*_m[1][0] + q[2][2]*_m[2][0];
  mout[2][1] = q[0][2]*_m[0][1] + q[1][2]*_m[1][1] + q[2][2]*_m[2][1];
  mout[2][2] = q[0][2]*_m[0][2] + q[1][2]*_m[1][2] + q[2][2]*_m[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute transformation from relative to absolute
 *        reference frame Q^t M Q
 *
 * \param[in]     m      symetric matrix of 3x3 real values
 * \param[in]     q      transformation matrix of 3x3 real values
 * \param[out]    mout   Q^t M Q
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_transform_r_to_a(const cs_real_t  m[6],
                                const cs_real_t  q[3][3],
                                cs_real_t        mout[6])
{
  /* _m = M.Q */
  cs_real_33_t _m;
  _m[0][0] = m[0]*q[0][0] + m[3]*q[1][0] + m[5]*q[2][0];
  _m[0][1] = m[0]*q[0][1] + m[3]*q[1][1] + m[5]*q[2][1];
  _m[0][2] = m[0]*q[0][2] + m[3]*q[1][2] + m[5]*q[2][2];

  _m[1][0] = m[3]*q[0][0] + m[1]*q[1][0] + m[4]*q[2][0];
  _m[1][1] = m[3]*q[0][1] + m[1]*q[1][1] + m[4]*q[2][1];
  _m[1][2] = m[3]*q[0][2] + m[1]*q[1][2] + m[4]*q[2][2];

  _m[2][0] = m[5]*q[0][0] + m[4]*q[1][0] + m[2]*q[2][0];
  _m[2][1] = m[5]*q[0][1] + m[4]*q[1][1] + m[2]*q[2][1];
  _m[2][2] = m[5]*q[0][2] + m[4]*q[1][2] + m[2]*q[2][2];

  /* mout = Q^t _m */
  mout[0] = q[0][0]*_m[0][0] + q[1][0]*_m[1][0] + q[2][0]*_m[2][0];
  mout[1] = q[0][1]*_m[0][1] + q[1][1]*_m[1][1] + q[2][1]*_m[2][1];
  mout[2] = q[0][2]*_m[0][2] + q[1][2]*_m[1][2] + q[2][2]*_m[2][2];

  mout[3] = q[0][0]*_m[0][1] + q[1][0]*_m[1][1] + q[2][0]*_m[2][1];
  mout[4] = q[0][1]*_m[0][2] + q[1][1]*_m[1][2] + q[2][1]*_m[2][2];
  mout[5] = q[0][0]*_m[0][2] + q[1][0]*_m[1][2] + q[2][0]*_m[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute transformation from absolute to relative
 *        reference frame Q M Q^t
 *
 * \param[in]     m      matrix of 3x3 real values
 * \param[in]     q      transformation matrix of 3x3 real values
 * \param[out]    mout   Q M Q^t
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_transform_a_to_r(const cs_real_t  m[3][3],
                            const cs_real_t  q[3][3],
                            cs_real_t        mout[3][3])
{
  /* _m = M.Q^t */
  cs_real_33_t _m;
  _m[0][0] = m[0][0]*q[0][0] + m[0][1]*q[0][1] + m[0][2]*q[0][2];
  _m[0][1] = m[0][0]*q[1][0] + m[0][1]*q[1][1] + m[0][2]*q[1][2];
  _m[0][2] = m[0][0]*q[2][0] + m[0][1]*q[2][1] + m[0][2]*q[2][2];

  _m[1][0] = m[1][0]*q[0][0] + m[1][1]*q[0][1] + m[1][2]*q[0][2];
  _m[1][1] = m[1][0]*q[1][0] + m[1][1]*q[1][1] + m[1][2]*q[1][2];
  _m[1][2] = m[1][0]*q[2][0] + m[1][1]*q[2][1] + m[1][2]*q[2][2];

  _m[2][0] = m[2][0]*q[0][0] + m[2][1]*q[0][1] + m[2][2]*q[0][2];
  _m[2][1] = m[2][0]*q[1][0] + m[2][1]*q[1][1] + m[2][2]*q[1][2];
  _m[2][2] = m[2][0]*q[2][0] + m[2][1]*q[2][1] + m[2][2]*q[2][2];

  /* mout = Q _m */
  mout[0][0] = q[0][0]*_m[0][0] + q[0][1]*_m[1][0] + q[0][2]*_m[2][0];
  mout[0][1] = q[0][0]*_m[0][1] + q[0][1]*_m[1][1] + q[0][2]*_m[2][1];
  mout[0][2] = q[0][0]*_m[0][2] + q[0][1]*_m[1][2] + q[0][2]*_m[2][2];

  mout[1][0] = q[1][0]*_m[0][0] + q[1][1]*_m[1][0] + q[1][2]*_m[2][0];
  mout[1][1] = q[1][0]*_m[0][1] + q[1][1]*_m[1][1] + q[1][2]*_m[2][1];
  mout[1][2] = q[1][0]*_m[0][2] + q[1][1]*_m[1][2] + q[1][2]*_m[2][2];

  mout[2][0] = q[2][0]*_m[0][0] + q[2][1]*_m[1][0] + q[2][2]*_m[2][0];
  mout[2][1] = q[2][0]*_m[0][1] + q[2][1]*_m[1][1] + q[2][2]*_m[2][1];
  mout[2][2] = q[2][0]*_m[0][2] + q[2][1]*_m[1][2] + q[2][2]*_m[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute transformation from absolute to relative
 *        reference frame Q M Q^t
 *
 * \param[in]     m      symetric matrix of 3x3 real values
 * \param[in]     q      transformation matrix of 3x3 real values
 * \param[out]    mout   Q M Q^t
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_transform_a_to_r(const cs_real_t  m[6],
                                const cs_real_t  q[3][3],
                                cs_real_t        mout[6])
{
  /* _m = M.Q^t */
  cs_real_33_t _m;
  _m[0][0] = m[0]*q[0][0] + m[3]*q[0][1] + m[5]*q[0][2];
  _m[0][1] = m[0]*q[1][0] + m[3]*q[1][1] + m[5]*q[1][2];
  _m[0][2] = m[0]*q[2][0] + m[3]*q[2][1] + m[5]*q[2][2];

  _m[1][0] = m[3]*q[0][0] + m[1]*q[0][1] + m[4]*q[0][2];
  _m[1][1] = m[3]*q[1][0] + m[1]*q[1][1] + m[4]*q[1][2];
  _m[1][2] = m[3]*q[2][0] + m[1]*q[2][1] + m[4]*q[2][2];

  _m[2][0] = m[5]*q[0][0] + m[4]*q[0][1] + m[2]*q[0][2];
  _m[2][1] = m[5]*q[1][0] + m[4]*q[1][1] + m[2]*q[1][2];
  _m[2][2] = m[5]*q[2][0] + m[4]*q[2][1] + m[2]*q[2][2];

  /* mout = Q _m */
  mout[0] = q[0][0]*_m[0][0] + q[0][1]*_m[1][0] + q[0][2]*_m[2][0];
  mout[1] = q[1][0]*_m[0][1] + q[1][1]*_m[1][1] + q[1][2]*_m[2][1];
  mout[2] = q[2][0]*_m[0][2] + q[2][1]*_m[1][2] + q[2][2]*_m[2][2];


  mout[3] = q[0][0]*_m[0][1] + q[0][1]*_m[1][1] + q[0][2]*_m[2][1];
  mout[4] = q[1][0]*_m[0][2] + q[1][1]*_m[1][2] + q[1][2]*_m[2][2];
  mout[5] = q[0][0]*_m[0][2] + q[0][1]*_m[1][2] + q[0][2]*_m[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Extract from the given matrix its symmetric
 *        and anti-symmetric part
 *
 * \param[in]     m        matrix of 3x3 real values
 * \param[out]    m_sym    matrix of 3x3 real values (symmetric part)
 * \param[out]    m_ant    matrix of 3x3 real values (anti-symmetric part)
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_extract_sym_ant(const cs_real_t  m[3][3],
                           cs_real_t        m_sym[3][3],
                           cs_real_t        m_ant[3][3])
{
  /* sym = 0.5 (m + m_transpose) */
  m_sym[0][0] = 0.5 * (m[0][0] + m[0][0]);
  m_sym[0][1] = 0.5 * (m[0][1] + m[1][0]);
  m_sym[0][2] = 0.5 * (m[0][2] + m[2][0]);
  m_sym[1][0] = 0.5 * (m[1][0] + m[0][1]);
  m_sym[1][1] = 0.5 * (m[1][1] + m[1][1]);
  m_sym[1][2] = 0.5 * (m[1][2] + m[2][1]);
  m_sym[2][0] = 0.5 * (m[2][0] + m[0][2]);
  m_sym[2][1] = 0.5 * (m[2][1] + m[1][2]);
  m_sym[2][2] = 0.5 * (m[2][2] + m[2][2]);

  /* ant = 0.5 (m - m_transpose) */
  m_ant[0][0] = 0.5 * (m[0][0] - m[0][0]);
  m_ant[0][1] = 0.5 * (m[0][1] - m[1][0]);
  m_ant[0][2] = 0.5 * (m[0][2] - m[2][0]);
  m_ant[1][0] = 0.5 * (m[1][0] - m[0][1]);
  m_ant[1][1] = 0.5 * (m[1][1] - m[1][1]);
  m_ant[1][2] = 0.5 * (m[1][2] - m[2][1]);
  m_ant[2][0] = 0.5 * (m[2][0] - m[0][2]);
  m_ant[2][1] = 0.5 * (m[2][1] - m[1][2]);
  m_ant[2][2] = 0.5 * (m[2][2] - m[2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the second main invariant of the symmetric part
 *         of a 3x3 tensor
 *
 * \param[in]     m        matrix of 3x3 real values
 *
 * \return trace(S.S) = Sij Sij, where S is the symmetric part of m
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline cs_real_t
cs_math_33_main_invariant_2(const cs_real_t  m[3][3])
{
  /* sym = 0.5 (m + m_transpose) */
  return cs_math_pow2(m[0][0])
       + cs_math_pow2(m[1][1])
       + cs_math_pow2(m[2][2])
    + 0.5 * ( cs_math_pow2(m[0][1] + m[1][0])
            + cs_math_pow2(m[0][2] + m[2][0])
            + cs_math_pow2(m[1][2] + m[2][1]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the product of two 3x3 real matrices to a matrix.
 *
 * \param[in]     m1            matrix of 3x3 real values
 * \param[in]     m2            matrix of 3x3 real values
 * \param[out]    mout          matrix of 3x3 real values
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_33_product_add(const cs_real_t             m1[3][3],
                       const cs_real_t             m2[3][3],
                       cs_real_t        (*restrict mout)[3])
{
  mout[0][0] += m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
  mout[0][1] += m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1];
  mout[0][2] += m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2];

  mout[1][0] += m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
  mout[1][1] += m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1];
  mout[1][2] += m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2];

  mout[2][0] += m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
  mout[2][1] += m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1];
  mout[2][2] += m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of two symmetric matrices.
 *
 * Warning: this is valid if and only if s1 and s2 commute (otherwise sout is
 *          not symmetric).
 *
 * \remark Symmetric matrix coefficients are stored as follows:
 *         (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s1            symmetric matrix
 * \param[in]     s2            symmetric matrix
 * \param[out]    sout          sout = s1 * s2
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_product(const cs_real_t      s1[6],
                       const cs_real_t      s2[6],
                       cs_real_t *restrict  sout)
{
  /* S11 */
  sout[0] = s1[0]*s2[0] + s1[3]*s2[3] + s1[5]*s2[5];
  /* S22 */
  sout[1] = s1[3]*s2[3] + s1[1]*s2[1] + s1[4]*s2[4];
  /* S33 */
  sout[2] = s1[5]*s2[5] + s1[4]*s2[4] + s1[2]*s2[2];
  /* S12 = S21 */
  sout[3] = s1[0]*s2[3] + s1[3]*s2[1] + s1[5]*s2[4];
  /* S23 = S32 */
  sout[4] = s1[3]*s2[5] + s1[1]*s2[4] + s1[4]*s2[2];
  /* S13 = S31 */
  sout[5] = s1[0]*s2[5] + s1[3]*s2[4] + s1[5]*s2[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a 6x6 matrix A, equivalent to a 3x3 matrix s, such as:
 *        A*R_6 = R*s^t + s*R
 *
 * \param[in]     s            3x3 matrix
 * \param[out]    sout         6x6 matrix
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_reduce_sym_prod_33_to_66(const cs_real_t       s[3][3],
                                 cs_real_t  (*restrict sout)[6])
{
  const int t2v[3][3] = {{0, 3, 5},
                         {3, 1, 4},
                         {5, 4, 2}};

  const int iv2t[6] = {0, 1, 2, 0, 1, 0};
  const int jv2t[6] = {0, 1, 2, 1, 2, 2};

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++)
      sout[i][j] = 0;
  }

  /* Consider : W = s*R + R*s^t .
   *            W_ij = Sum_{k<3} [s_ik*r_jk + s_jk*r_ik]
   * We look for A_(ij,pq) such as A*R = W
   *
   * so
   *   A_(ij,jk) takes s_ik
   * and
   *   A_(ij,ik) takes s_jk
   */
  for (int ij = 0; ij < 6; ij++) {
    int i = iv2t[ij];
    int j = jv2t[ij];
    for (int k = 0; k < 3; k++) {
      int ik = t2v[i][k];
      int jk = t2v[j][k];

      sout[ij][ik] += s[j][k];
      sout[ij][jk] += s[i][k];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of three symmetric matrices.
 *
 * \remark Symmetric matrix coefficients are stored as follows:
 *         (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s1            symmetric matrix
 * \param[in]     s2            symmetric matrix
 * \param[in]     s3            symmetric matrix
 * \param[out]    sout          sout = s1 * s2 * s3
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_math_sym_33_double_product(const cs_real_t       s1[6],
                              const cs_real_t       s2[6],
                              const cs_real_t       s3[6],
                              cs_real_t  (*restrict sout)[3])
{
  cs_real_t _sout[3][3];

  /* S11 */
  _sout[0][0] = s1[0]*s2[0] + s1[3]*s2[3] + s1[5]*s2[5];
  /* S22 */
  _sout[1][1] = s1[3]*s2[3] + s1[1]*s2[1] + s1[4]*s2[4];
  /* S33 */
  _sout[2][2] = s1[5]*s2[5] + s1[4]*s2[4] + s1[2]*s2[2];
  /* S12  */
  _sout[0][1] = s1[0]*s2[3] + s1[3]*s2[1] + s1[5]*s2[4];
  /* S21  */
  _sout[1][0] = s2[0]*s1[3] + s2[3]*s1[1] + s2[5]*s1[4];
  /* S23  */
  _sout[1][2] = s1[3]*s2[5] + s1[1]*s2[4] + s1[4]*s2[2];
  /* S32  */
  _sout[2][1] = s2[3]*s1[5] + s2[1]*s1[4] + s2[4]*s1[2];
  /* S13  */
  _sout[0][2] = s1[0]*s2[5] + s1[3]*s2[4] + s1[5]*s2[2];
  /* S31  */
  _sout[2][0] = s2[0]*s1[5] + s2[3]*s1[4] + s2[5]*s1[2];

  sout[0][0] = _sout[0][0]*s3[0] + _sout[0][1]*s3[3] + _sout[0][2]*s3[5];
  /* S22 */
  sout[1][1] = _sout[1][0]*s3[3] + _sout[1][1]*s3[1] + _sout[1][2]*s3[4];
  /* S33 */
  sout[2][2] = _sout[2][0]*s3[5] + _sout[2][1]*s3[4] + _sout[2][2]*s3[2];
  /* S12  */
  sout[0][1] = _sout[0][0]*s3[3] + _sout[0][1]*s3[1] + _sout[0][2]*s3[4];
  /* S21  */
  sout[1][0] = s3[0]*_sout[1][0] + s3[3]*_sout[1][1] + s3[5]*_sout[1][2];
  /* S23  */
  sout[1][2] = _sout[1][0]*s3[5] + _sout[1][1]*s3[4] + _sout[1][2]*s3[2];
  /* S32  */
  sout[2][1] = s3[3]*_sout[2][0] + s3[1]*_sout[2][1] + s3[4]*_sout[2][2];
  /* S13  */
  sout[0][2] = _sout[0][0]*s3[5] + _sout[0][1]*s3[4] + _sout[0][2]*s3[2];
  /* S31  */
  sout[2][0] = s3[0]*_sout[2][0] + s3[3]*_sout[2][1] + s3[5]*_sout[2][2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_nvec3_t structure from a cs_real_3_t
 *
 * \param[in]  v     vector of size 3
 * \param[out] qv    pointer to a cs_nvec3_t structure
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static inline void
cs_nvec3(const cs_real_3_t    v,
         cs_nvec3_t          *qv)
{
  cs_real_t  magnitude = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  qv->meas = magnitude;
  if (fabs(magnitude) > cs_math_zero_threshold) {

    const cs_real_t  inv = 1/magnitude;
    qv->unitv[0] = inv * v[0];
    qv->unitv[1] = inv * v[1];
    qv->unitv[2] = inv * v[2];

  }
  else
    qv->unitv[0] = qv->unitv[1] = qv->unitv[2] = 0;
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief  Compute the length (Euclidean norm) between two points xa and xb in
 *         a Cartesian coordinate system of dimension 3
 *
 * \param[in]   xa       coordinate of the first extremity
 * \param[in]   xb       coordinate of the second extremity
 * \param[out]  len      pointer to the length of the vector va -> vb
 * \param[out]  unitv    unitary vector along xa -> xb
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE void
cs_math_3_length_unitv(const cs_real_t    xa[3],
                       const cs_real_t    xb[3],
                       cs_real_t         *len,
                       cs_real_3_t        unitv);

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
                     cs_real_t        eig_vals[3]);

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
                 cs_real_t          *eig_max);

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
                const cs_real_t  xf[3]);

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
               const cs_real_t   xc[3]);

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
cs_math_33_eig_val_vec(const cs_real_t  m_in[3][3],
                       const cs_real_t  tol_err,
                       cs_real_t        eig_val[3],
                       cs_real_t        eig_vec[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute LU factorization of an array of dense matrices
 *        of identical size.
 *
 * \param[in]   n         matrix size
 * \param[in]   a         matrix (size n.n)
 * \param[out]  a_lu      LU factorization of matrix (size n.n)
 */
/*----------------------------------------------------------------------------*/

void
cs_math_fact_lu(const int         n,
                const cs_real_t  *a,
                cs_real_t        *a_lu);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute forward and backward to solve an LU P*P system.
 *
 * \param[in]   a_lu   matrix LU factorization
 * \param[in]   n      matrix size
 * \param[out]  x      solution
 * \param[out]  b      right hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_math_fw_and_bw_lu(const cs_real_t  a_lu[],
                     const int        n,
                     cs_real_t        x[],
                     const cs_real_t  b[]);

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
cs_math_sym_44_factor_ldlt(cs_real_t  ldlt[10]);

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
                                  const cs_real_t  rhs[4]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATH_H__ */
