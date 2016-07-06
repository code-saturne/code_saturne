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

static cs_real_t  _machine_epsilon = 1.11e-16;

/* Numerical constants */

const cs_real_t cs_math_zero_threshold = FLT_MIN;
const cs_real_t cs_math_onethird = 1./3.;
const cs_real_t cs_math_onesix = 1./6.;
const cs_real_t cs_math_onetwelve = 1./12.;

/*! epsilon \f$ 10^{-12}\f$ */
const cs_real_t cs_math_epzero = 1e-12;

/*! infinite \f$ 10^{+30}\f$ */
const cs_real_t cs_math_infinite_r = 1.e30;

/*! big value \f$ 10^{+12}\f$ */
const cs_real_t cs_math_big_r = 1.e12;

/*! \f$ \pi \f$ value with 20 digits */
const cs_real_t cs_math_pi = 3.14159265358979323846;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_inv_cramer
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_inverse, SYMMETRIC_MATRIX_INVERSE)
(
  const cs_real_6_t s,
  cs_real_6_t       sout
)
{
  cs_math_sym_33_inv_cramer(s,
                            sout);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_math_sym_33_product
 *----------------------------------------------------------------------------*/

void CS_PROCF (symmetric_matrix_product, SYMMETRIC_MATRIX_PRODUCT)
(
  const cs_real_6_t s1,
  const cs_real_6_t s2,
  cs_real_6_t       sout
)
{
  cs_math_sym_33_product(s1,
                         s2,
                         sout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to the machine precision
 */
/*----------------------------------------------------------------------------*/

void
cs_math_set_machine_epsilon(void)
{
  double  eps = 5e-16;
  double  y = 1.0 + eps;

  while (y > 1.0) {
    eps /= 2.0;
    y = 1.0 + eps;
  }
  eps *= 2.0;

  _machine_epsilon = eps;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the value related to the machine precision
 */
/*----------------------------------------------------------------------------*/

double
cs_math_get_machine_epsilon(void)
{
  return _machine_epsilon;
}

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
cs_math_33_eigen(const cs_real_t     m[3][3],
                 cs_real_t          *eig_ratio,
                 cs_real_t          *eig_max)
{
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

  if (p1 > 0.0) { // m is not diagonal

    cs_real_t  theta;
    cs_real_t  n[3][3];
    cs_real_t  tr = cs_math_onethird*(m[0][0] + m[1][1] + m[2][2]);

    e1 = m[0][0] - tr, e2 = m[1][1] - tr, e3 = m[2][2] - tr;
    cs_real_t  p2 = e1*e1 + e2*e2 + e3*e3 + 2*p1;

    assert(p2 > 0);
    cs_real_t  p = sqrt(p2*cs_math_onesix);
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
    cs_real_t  r = 0.5 * cs_math_33_determinant((const cs_real_t (*)[3])n);

    if (r <= -1)
      theta = cs_math_onethird*cs_math_pi;
    else if (r >= 1)
      theta = 0.;
    else
      theta = cs_math_onethird*acos(r);

    // eigenvalues computed should satisfy e1 < e2 < e3
    e3 = tr + 2*p*cos(theta);
    e1 = tr + 2*p*cos(theta + 2*cs_math_pi*cs_math_onethird);
    e2 = 3*tr - e1 -e3; // since tr(m) = e1 + e2 + e3

  }
  else { // m is diagonal

    e1 = m[0][0], e2 = m[1][1], e3 = m[2][2];
    if (e3 < e2) e = e3, e3 = e2, e2 = e;
    if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
    else {
      if (e2 < e1) e = e2, e2 = e1, e1 = e;
    }

  } /* diagonal or not */

  /* Return values */
  assert(fabs(e1) > 0);
  *eig_ratio = e3/e1;
  *eig_max = e3;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the length (euclidien norm) between two points xa and xb in
 *         a cartesian coordinate system of dimension 3
 *
 * \param[in]   xa       coordinate of the first extremity
 * \param[in]   xb       coordinate of the second extremity
 * \param[out]  len      pointer to the length of the vector va -> vb
 * \param[out]  unitv    unitary vector anlong va -> vb
 */
/*----------------------------------------------------------------------------*/

inline void
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

inline double
cs_math_surftri(const cs_real_t  xv[3],
                const cs_real_t  xe[3],
                const cs_real_t  xf[3])
{
  cs_real_3_t  u, v, cp;

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
  double  lev, lef, lec;
  cs_real_3_t  uev, uef, uec, ucp;

  cs_math_3_length_unitv(xe, xv, &lev, uev);
  cs_math_3_length_unitv(xe, xf, &lef, uef);
  cs_math_3_length_unitv(xe, xc, &lec, uec);
  cs_math_3_cross_product(uev, uef, ucp);

  return  cs_math_onesix *lev*lef*lec* fabs(cs_math_3_dot_product(ucp, uec));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
