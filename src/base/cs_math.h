#ifndef __CS_MATH_H__
#define __CS_MATH_H__

/*============================================================================
 * Mathematical base functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_inv_cramer
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_inverse, SYMMETRIC_MATRIX_INVERSE)
(
  const cs_real_6_t s,
  cs_real_6_t       sout
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_product
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_product, SYMMETRIC_MATRIX_PRODUCT)
(
 const cs_real_6_t s1,
 const cs_real_6_t s2,
 cs_real_6_t       sout
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of a matrix of 3x3 real values by a vector of 3
 * real values.
 *
 * \param[in]     m             matrix of 3x3 real values
 * \param[in]     v             vector of 3 real values
 * \param[out]    mv            vector of 3 real values
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_math_33_3_product(const cs_real_33_t m,
                     const cs_real_3_t  v,
                     cs_real_3_t mv)
{
  for (int ii = 0; ii < 3; ii++)
    mv[ii] = m[ii][0] * v[0] + m[ii][1] * v[1] + m[ii][2] * v[2];
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
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_math_sym_33_3_product(const cs_real_6_t  m,
                         const cs_real_3_t  v,
                         cs_real_3_t        mv)
{
  mv[0] = m[0] * v[0] + m[3] * v[1] + m[5] * v[2];
  mv[1] = m[3] * v[0] + m[1] * v[1] + m[4] * v[2];
  mv[2] = m[5] * v[0] + m[4] * v[1] + m[2] * v[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dot product of two vectors of 3 real values.
 *
 * \param[in]     u             vector of 3 real values
 * \param[in]     v             vector of 3 real values
 *
 * \return the resulting dot product u.v.
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_math_3_dot_product(const cs_real_3_t u,
                      const cs_real_3_t v)
{
  cs_real_t uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

  return uv;
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

static inline cs_real_t
cs_math_3_square_norm(const cs_real_3_t v)
{
  cs_real_t v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

  return v2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of a symmetric matrix using Cramer's rule.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s             symmetric matrix
 * \param[out]    sout          sout = s1 * s2
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_math_sym_33_inv_cramer(const cs_real_6_t s,
                          cs_real_6_t       sout)
{
  double detinv;

  sout[0] = s[1]*s[2] - s[4]*s[4];
  sout[1] = s[0]*s[2] - s[5]*s[5];
  sout[2] = s[0]*s[1] - s[3]*s[3];
  sout[3] = s[4]*s[5] - s[3]*s[2];
  sout[4] = s[3]*s[5] - s[0]*s[4];
  sout[5] = s[3]*s[4] - s[1]*s[5];

  detinv = 1. / (s[0]*sout[0] + s[3]*sout[3] + s[5]*sout[5]);

  sout[0] = sout[0] * detinv;
  sout[1] = sout[1] * detinv;
  sout[2] = sout[2] * detinv;
  sout[3] = sout[3] * detinv;
  sout[4] = sout[4] * detinv;
  sout[5] = sout[5] * detinv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the product of two symmetric matrices.
 * NB: Symmetric matrix are stored as follows (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s1            symmetric matrix
 * \param[in]     s2            symmetric matrix
 * \param[out]    sout          sout = s1 * s2
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_math_sym_33_product(const cs_real_6_t s1,
                       const cs_real_6_t s2,
                       cs_real_6_t       sout)
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

END_C_DECLS

#endif /* __CS_MATH_H__ */
