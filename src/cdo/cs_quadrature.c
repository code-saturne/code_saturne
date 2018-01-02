/*============================================================================
 * Routines to handle quadrature rules
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"

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
 * \param[in, out]  w        weight (same value for all points)
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_4pts(const cs_real_3_t  xv,
                       const cs_real_3_t  xe,
                       const cs_real_3_t  xf,
                       const cs_real_3_t  xc,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double            *w)
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
  *w = 0.25*vol;
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
  for (int k = 0; k < 4; k++)
    weights[k] = wv2;
  weights[4] = wv1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return th name associated to a type of quadrature
 *
 * \param[in]     type     cs_quadra_type_t
 *
 * \return the name associated to a given type of quadrature
 */
/*----------------------------------------------------------------------------*/

const char *
cs_quadrature_get_type_name(const cs_quadra_type_t  type)
{
  return cs_quadrature_type_name[type];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
