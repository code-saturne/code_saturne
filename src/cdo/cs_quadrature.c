/*============================================================================
 * Routines to handle quadrature rules
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_quadrature.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

static const double  _quad_25ov48 = 25./48.;
static const double  _quad_9ov16 = 9./16.;
static const double  _quad_over18 = 1./18.;
static const double  _quad_over6 = 1./6.;
static const double  _quad_over3 = 1./3.;
static const double  _quad_9ov40 = 9./40.;
static const double  _quad_31ov240 = 31./240.;
static const double  _tetr_quad15w3 = 10. / 189.;
static const double  _tetr_quad15w4 = 16. / 135.;

/* Constant quadrature weights */
static double  _edge_quad2c1;
static double  _edge_quad2c2;
static double  _edge_quad3c1;
static double  _edge_quad3c2;
static double  _tria_quad7c1;
static double  _tria_quad7c2;
static double  _tria_quad7c3;
static double  _tria_quad7c4;
static double  _tetr_quad4c1;
static double  _tetr_quad4c2;
static double  _tetr_quad15g1;
static double  _tetr_quad15g11;
static double  _tetr_quad15g2;
static double  _tetr_quad15g21;
static double  _tetr_quad15g3;
static double  _tetr_quad15g31;
static double  _tetr_quad15w1;
static double  _tetr_quad15w2;

static const char
cs_quadrature_type_name[CS_QUADRATURE_N_TYPES][CS_BASE_STRING_LEN] =
  { N_("none"),
    N_("barycentric"),
    N_("barycentric on a tetrahedral submesh"),
    N_("higher (single weight)"),
    N_("highest") };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute constant weights for all quadratures used
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_setup(void)
{
  /* Quadrature on edge with 2 points */
  _edge_quad2c1 = 0.5*(1 + _quad_over3*sqrt(3.));
  _edge_quad2c2 = 0.5*(1 - _quad_over3*sqrt(3.));

  /* Quadrature on edge with 3 points */
  _edge_quad3c1 = 0.5*(1 + sqrt(0.6));
  _edge_quad3c2 = 0.5*(1 - sqrt(0.6));

  /* Quadrature on triangles with 7 points */
  _tria_quad7c1 = (6. - sqrt(15.))/21.;
  _tria_quad7c2 = (6. + sqrt(15.))/21.;
  _tria_quad7c3 = _quad_31ov240 - sqrt(15.)/1200.;
  _tria_quad7c4 = _quad_31ov240 + sqrt(15.)/1200.;

  /* Quadrature on tetrahedron with 4 points */
  _tetr_quad4c1 = 0.05*(5. - sqrt(5));
  _tetr_quad4c2 = 1. -3.*_tetr_quad4c1;

  /* Quadrature on tetrahedron with 15 points */
  _tetr_quad15g1  = (7. - sqrt(15.) ) / 34.;
  _tetr_quad15g11 = 1. - 3. * _tetr_quad15g1;
  _tetr_quad15g2  = 7./17. - _tetr_quad15g1;
  _tetr_quad15g21 = 1. - 3. * _tetr_quad15g2;
  _tetr_quad15g3  = (5. - sqrt(15.)) / 20.;
  _tetr_quad15g31 = (5. + sqrt(15.)) / 20.;
  _tetr_quad15w1  = (2665. + 14. * sqrt(15.) ) / 37800.;
  _tetr_quad15w2  = (2665. - 14. * sqrt(15.) ) / 37800.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for an edge from v1 -> v2 (2 points)
 *          Exact for polynomial function up to order 3
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      len      length of edge [v1, v2]
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weight (same weight for the two points)
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_edge_2pts(const cs_real_3_t  v1,
                        const cs_real_3_t  v2,
                        double             len,
                        cs_real_3_t        gpts[],
                        double            *w)
{
  int  k;

  /* Compute quadrature points */
  for (k = 0; k < 3; k++) {
    gpts[0][k] = v1[k]*_edge_quad2c2 + v2[k]*_edge_quad2c1;
    gpts[1][k] = v1[k]*_edge_quad2c1 + v2[k]*_edge_quad2c2;
  }

  /* Compute weights */
  *w= 0.5*len;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for an edge from v1 -> v2 (3 points)
 *          Exact for polynomial function up to order 5
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      len      length of edge [v1, v2]
 * \param[in ,out] gpts     gauss points
 * \param[in, out] w        weights
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_edge_3pts(const cs_real_3_t  v1,
                        const cs_real_3_t  v2,
                        double             len,
                        cs_real_3_t        gpts[],
                        double             w[])
{
  int  k;

  const double  b = len * _quad_over18;

  /* Compute quadrature points */
  for (k = 0; k < 3; k++) {
    gpts[0][k] = 0.5*( v1[k] + v2[k] );
    gpts[1][k] = v1[k]*_edge_quad3c1 + v2[k]*_edge_quad3c2;
    gpts[2][k] = v1[k]*_edge_quad3c2 + v2[k]*_edge_quad3c1;
  }

  /* Compute weights */
  w[0]= 8*b, w[1] = w[2] = 5*b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for a triangle (3 points)
 *          Exact for polynomial function up to order 2
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      v3       third vertex
 * \param[in]      area     area of triangle {v1, v2, v3}
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weight (same weight for the three points)
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tria_3pts(const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        const cs_real_3_t   v3,
                        double              area,
                        cs_real_3_t         gpts[],
                        double             *w)
{
  int  k;

  /* Compute quadrature points */
  for (k = 0; k < 3; k++) {
    gpts[0][k] = 0.5*( v1[k] + v2[k] );
    gpts[1][k] = 0.5*( v1[k] + v3[k] );
    gpts[2][k] = 0.5*( v2[k] + v3[k] );
  }

  /* Compute weight */
  *w = _quad_over3 * area;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for a triangle (4 points)
 *          Exact for polynomial function up to order 3
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      v3       third vertex
 * \param[in]      area     area of triangle {v1, v2, v3}
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weights
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tria_4pts(const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        const cs_real_3_t   v3,
                        double              area,
                        cs_real_3_t         gpts[],
                        double              w[])
{
  int  k;

  /* Compute quadrature points */
  for (k = 0; k < 3; k++) {
    gpts[0][k] = _quad_over3*( v1[k] + v2[k] + v3[k] );
    gpts[1][k] = 0.2*( v1[k] + v2[k] ) + 0.6*v3[k];
    gpts[2][k] = 0.2*( v1[k] + v3[k] ) + 0.6*v2[k];
    gpts[3][k] = 0.2*( v2[k] + v3[k] ) + 0.6*v1[k];
  }

  /* Compute weight */
  w[0] = -_quad_9ov16 * area;
  w[1] = w[2] = w[3] = _quad_25ov48 * area;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for a triangle (7 points)
 *          Exact for polynomial function up to order 5
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      v3       third vertex
 * \param[in]      area     area of triangle {v1, v2, v3}
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weights
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tria_7pts(const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        const cs_real_3_t   v3,
                        double              area,
                        cs_real_3_t         gpts[],
                        double              w[])
{
  int  k;

  /* Compute quadrature points */
  for (k = 0; k < 3; k++) {
    gpts[0][k] = _quad_over3*( v1[k] + v2[k] + v3[k] );
    gpts[1][k] = _tria_quad7c1* (v1[k] + v2[k]) + (1-2*_tria_quad7c1)*v3[k];
    gpts[2][k] = _tria_quad7c1* (v3[k] + v1[k]) + (1-2*_tria_quad7c1)*v2[k];
    gpts[3][k] = _tria_quad7c1* (v2[k] + v3[k]) + (1-2*_tria_quad7c1)*v1[k];
    gpts[4][k] = _tria_quad7c2* (v1[k] + v2[k]) + (1-2*_tria_quad7c2)*v3[k];
    gpts[5][k] = _tria_quad7c2* (v3[k] + v1[k]) + (1-2*_tria_quad7c2)*v2[k];
    gpts[6][k] = _tria_quad7c2* (v2[k] + v3[k]) + (1-2*_tria_quad7c2)*v1[k];
  }

  /* Compute weights */
  w[0] = _quad_9ov40 * area;
  w[1] = w[2] = w[3] = _tria_quad7c3 * area;
  w[4] = w[5] = w[6] = _tria_quad7c4 * area;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the quadrature in a tetrehedra. Exact for 2nd order
 *         polynomials (order 3).
 *
 * \param[in]       xv       first vertex
 * \param[in]       xe       second vertex
 * \param[in]       xf       third vertex
 * \param[in]       xc       fourth vertex
 * \param[in]       vol      volume of tetrahedron {xv, xe, xf, xc}
 * \param[in, out]  gpts     4 Gauss points (size = 3*4)
 * \param[in, out]  weights  weight (same value for all points)
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_4pts(const cs_real_3_t  xv,
                       const cs_real_3_t  xe,
                       const cs_real_3_t  xf,
                       const cs_real_3_t  xc,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double             weights[])
{
  /* Compute Gauss points */
  for (int k = 0; k < 3; k++) {

    const double xve = xv[k] + xe[k], xfc = xf[k] + xc[k];

    gpts[0][k] = _tetr_quad4c1*(xf[k] + xve) + _tetr_quad4c2*xc[k];
    gpts[1][k] = _tetr_quad4c1*(xe[k] + xfc) + _tetr_quad4c2*xv[k];
    gpts[2][k] = _tetr_quad4c1*(xv[k] + xfc) + _tetr_quad4c2*xe[k];
    gpts[3][k] = _tetr_quad4c1*(xc[k] + xve) + _tetr_quad4c2*xf[k];

  }

  /* Compute weights (multiplicity 4) */
  weights[0] = weights[1] = weights[2] = weights[3] = 0.25*vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the quadrature in a tetrehedra. Exact for 3rd order
 *         polynomials (order 4).
 *
 * \param[in]       xv       first vertex
 * \param[in]       xe       second vertex
 * \param[in]       xf       third vertex
 * \param[in]       xc       fourth vertex
 * \param[in]       vol      volume of tetrahedron {xv, xe, xf, xc}
 * \param[in, out]  gpts     5 Gauss points (size = 3*5)
 * \param[in, out]  weights  5 weigths related to each Gauss point
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_5pts(const cs_real_3_t  xv,
                       const cs_real_3_t  xe,
                       const cs_real_3_t  xf,
                       const cs_real_3_t  xc,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double             weights[])
{
  const double  wv1 = -0.8*vol;  /* multiplicity 1 */
  const double  wv2 = 0.45*vol;  /* multiplicity 4 */

  /* compute Gauss points */
  for (int k = 0; k < 3; k++) {

    const double xve = xv[k] + xe[k], xfc = xf[k] + xc[k];

    gpts[0][k] = _quad_over6*(xve + xf[k]) + 0.5*xc[k];
    gpts[1][k] = _quad_over6*(xfc + xe[k]) + 0.5*xv[k];
    gpts[2][k] = _quad_over6*(xfc + xv[k]) + 0.5*xe[k];
    gpts[3][k] = _quad_over6*(xve + xc[k]) + 0.5*xf[k];
    gpts[4][k] = 0.25*(xve + xfc);
  }

  /* Compute weights */
  weights[0] = weights[1] = weights[2] = weights[3] = wv2;
  weights[4] = wv1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the quadrature in a tetrehedra. Exact for 5th order
 *         polynomials (order 6).
 *
 * \param[in]       xv       first vertex
 * \param[in]       xe       second vertex
 * \param[in]       xf       third vertex
 * \param[in]       xc       fourth vertex
 * \param[in]       vol      volume of tetrahedron {xv, xe, xf, xc}
 * \param[in, out]  gpts     15 Gauss points (size = 3*15)
 * \param[in, out]  weights  15 weigths related to each Gauss point
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_15pts(const cs_real_3_t   xv,
                        const cs_real_3_t   xe,
                        const cs_real_3_t   xf,
                        const cs_real_3_t   xc,
                        double              vol,
                        cs_real_3_t         gpts[],
                        double              weights[])
{
  const double w1 = vol * _tetr_quad15w1;
  const double w2 = vol * _tetr_quad15w2;
  const double w3 = vol * _tetr_quad15w3;

  for (short int i = 0; i < 3; ++i) {

    const double  xvxe = xv[i] + xe[i];
    const double  xvxf = xv[i] + xf[i];
    const double  xvxc = xv[i] + xc[i];
    const double  xexf = xe[i] + xf[i];
    const double  xexc = xe[i] + xc[i];
    const double  xfxc = xf[i] + xc[i];

    gpts[0][i]  = _tetr_quad15g1 * (xvxe + xf[i]) + _tetr_quad15g11 * xc[i];
    gpts[1][i]  = _tetr_quad15g1 * (xvxe + xc[i]) + _tetr_quad15g11 * xf[i];
    gpts[2][i]  = _tetr_quad15g1 * (xvxf + xc[i]) + _tetr_quad15g11 * xe[i];
    gpts[3][i]  = _tetr_quad15g1 * (xexf + xc[i]) + _tetr_quad15g11 * xv[i];

    gpts[4][i]  = _tetr_quad15g2 * (xvxe + xf[i]) + _tetr_quad15g21 * xc[i];
    gpts[5][i]  = _tetr_quad15g2 * (xvxe + xc[i]) + _tetr_quad15g21 * xf[i];
    gpts[6][i]  = _tetr_quad15g2 * (xvxf + xc[i]) + _tetr_quad15g21 * xe[i];
    gpts[7][i]  = _tetr_quad15g2 * (xexf + xc[i]) + _tetr_quad15g21 * xv[i];

    gpts[8][i]  = _tetr_quad15g3 * xvxe + _tetr_quad15g31 * xfxc;
    gpts[9][i]  = _tetr_quad15g3 * xvxc + _tetr_quad15g31 * xexf;
    gpts[10][i] = _tetr_quad15g3 * xvxf + _tetr_quad15g31 * xexc;
    gpts[11][i] = _tetr_quad15g3 * xexf + _tetr_quad15g31 * xvxc;
    gpts[12][i] = _tetr_quad15g3 * xfxc + _tetr_quad15g31 * xvxe;
    gpts[13][i] = _tetr_quad15g3 * xexc + _tetr_quad15g31 * xvxf;
    gpts[14][i] = 0.25* (xvxe + xfxc);

  }

  weights[0]  = weights[1] = weights[ 2] = weights[ 3] = w1;
  weights[4]  = weights[5] = weights[ 6] = weights[ 7] = w2;
  weights[8]  = weights[9] = weights[10] = weights[11] = weights[12] = w3;
  weights[13] = w3;
  weights[14] = vol * _tetr_quad15w4;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return th name associated to a type of quadrature
 *
 * \param[in]     type     cs_quadrature_type_t
 *
 * \return the name associated to a given type of quadrature
 */
/*----------------------------------------------------------------------------*/

const char *
cs_quadrature_get_type_name(const cs_quadrature_type_t  type)
{
  return cs_quadrature_type_name[type];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
