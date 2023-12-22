/*============================================================================
 * Functions to handle quadrature rules
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_flag.h"

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

/* Constant quadrature weights which are pre-computed */

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
 * \brief Compute constant weights for all quadratures used
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
 * \brief Return th name associated to a type of quadrature
 *
 * \param[in] type     cs_quadrature_type_t
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
                        double             w[])
{
  /* Compute quadrature points */

  for (int k = 0; k < 3; k++) {
    gpts[0][k] = v1[k]*_edge_quad2c2 + v2[k]*_edge_quad2c1;
    gpts[1][k] = v1[k]*_edge_quad2c1 + v2[k]*_edge_quad2c2;
  }

  /* Compute weights */

  w[0] = w[1] = 0.5*len;
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
  const double  b = len * _quad_over18;

  /* Compute quadrature points */

  for (int k = 0; k < 3; k++) {
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
                        double              w[])
{
  /* Compute quadrature points */

  for (int k = 0; k < 3; k++) {
    gpts[0][k] = 0.5*( v1[k] + v2[k] );
    gpts[1][k] = 0.5*( v1[k] + v3[k] );
    gpts[2][k] = 0.5*( v2[k] + v3[k] );
  }

  /* Compute weight */

  w[0] = w[1] = w[2] = _quad_over3 * area;
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
  /* Compute quadrature points */

  for (int k = 0; k < 3; k++) {
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
  /* Compute quadrature points */

  for (int k = 0; k < 3; k++) {
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
 * \param[in]       v1       first vertex
 * \param[in]       v2       second vertex
 * \param[in]       v3       third vertex
 * \param[in]       v4       fourth vertex
 * \param[in]       vol      volume of tetrahedron {v1, v2, v3, v4}
 * \param[in, out]  gpts     4 Gauss points (size = 3*4)
 * \param[in, out]  weights  weight (same value for all points)
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_4pts(const cs_real_3_t  v1,
                       const cs_real_3_t  v2,
                       const cs_real_3_t  v3,
                       const cs_real_3_t  v4,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double             weights[])
{
  /* Compute Gauss points */

  for (int k = 0; k < 3; k++) {

    const double v12 = v1[k] + v2[k], v34 = v3[k] + v4[k];

    gpts[0][k] = _tetr_quad4c1*(v3[k] + v12) + _tetr_quad4c2*v4[k];
    gpts[1][k] = _tetr_quad4c1*(v2[k] + v34) + _tetr_quad4c2*v1[k];
    gpts[2][k] = _tetr_quad4c1*(v1[k] + v34) + _tetr_quad4c2*v2[k];
    gpts[3][k] = _tetr_quad4c1*(v4[k] + v12) + _tetr_quad4c2*v3[k];

  }

  /* Compute weights (multiplicity 4) */

  weights[0] = weights[1] = weights[2] = weights[3] = 0.25*vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the quadrature in a tetrehedra. Exact for 3rd order
 *         polynomials (order 4).
 *
 * \param[in]       v1       first vertex
 * \param[in]       v2       second vertex
 * \param[in]       v3       third vertex
 * \param[in]       v4       fourth vertex
 * \param[in]       vol      volume of tetrahedron {v1, v2, v3, v4}
 * \param[in, out]  gpts     5 Gauss points (size = 3*5)
 * \param[in, out]  weights  5 weigths related to each Gauss point
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_5pts(const cs_real_3_t  v1,
                       const cs_real_3_t  v2,
                       const cs_real_3_t  v3,
                       const cs_real_3_t  v4,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double             weights[])
{
  const double  wv1 = -0.8*vol;  /* multiplicity 1 */
  const double  wv2 = 0.45*vol;  /* multiplicity 4 */

  /* compute Gauss points */

  for (int k = 0; k < 3; k++) {

    const double v12 = v1[k] + v2[k], v34 = v3[k] + v4[k];

    gpts[0][k] = _quad_over6*(v12 + v3[k]) + 0.5*v4[k];
    gpts[1][k] = _quad_over6*(v34 + v2[k]) + 0.5*v1[k];
    gpts[2][k] = _quad_over6*(v34 + v1[k]) + 0.5*v2[k];
    gpts[3][k] = _quad_over6*(v12 + v4[k]) + 0.5*v3[k];
    gpts[4][k] = 0.25*(v12 + v34);
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
 * \param[in]       v1       first vertex
 * \param[in]       v2       second vertex
 * \param[in]       v3       third vertex
 * \param[in]       v4       fourth vertex
 * \param[in]       vol      volume of tetrahedron {v1, v2, v3, v4}
 * \param[in, out]  gpts     15 Gauss points (size = 3*15)
 * \param[in, out]  weights  15 weigths related to each Gauss point
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_tet_15pts(const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        const cs_real_3_t   v3,
                        const cs_real_3_t   v4,
                        double              vol,
                        cs_real_3_t         gpts[],
                        double              weights[])
{
  const double w1 = vol * _tetr_quad15w1;
  const double w2 = vol * _tetr_quad15w2;
  const double w3 = vol * _tetr_quad15w3;

  for (short int i = 0; i < 3; ++i) {

    const double  v1v2 = v1[i] + v2[i];
    const double  v1v3 = v1[i] + v3[i];
    const double  v1v4 = v1[i] + v4[i];
    const double  v2v3 = v2[i] + v3[i];
    const double  v2v4 = v2[i] + v4[i];
    const double  v3v4 = v3[i] + v4[i];

    gpts[0][i]  = _tetr_quad15g1 * (v1v2 + v3[i]) + _tetr_quad15g11 * v4[i];
    gpts[1][i]  = _tetr_quad15g1 * (v1v2 + v4[i]) + _tetr_quad15g11 * v3[i];
    gpts[2][i]  = _tetr_quad15g1 * (v1v3 + v4[i]) + _tetr_quad15g11 * v2[i];
    gpts[3][i]  = _tetr_quad15g1 * (v2v3 + v4[i]) + _tetr_quad15g11 * v1[i];

    gpts[4][i]  = _tetr_quad15g2 * (v1v2 + v3[i]) + _tetr_quad15g21 * v4[i];
    gpts[5][i]  = _tetr_quad15g2 * (v1v2 + v4[i]) + _tetr_quad15g21 * v3[i];
    gpts[6][i]  = _tetr_quad15g2 * (v1v3 + v4[i]) + _tetr_quad15g21 * v2[i];
    gpts[7][i]  = _tetr_quad15g2 * (v2v3 + v4[i]) + _tetr_quad15g21 * v1[i];

    gpts[8][i]  = _tetr_quad15g3 * v1v2 + _tetr_quad15g31 * v3v4;
    gpts[9][i]  = _tetr_quad15g3 * v1v4 + _tetr_quad15g31 * v2v3;
    gpts[10][i] = _tetr_quad15g3 * v1v3 + _tetr_quad15g31 * v2v4;
    gpts[11][i] = _tetr_quad15g3 * v2v3 + _tetr_quad15g31 * v1v4;
    gpts[12][i] = _tetr_quad15g3 * v3v4 + _tetr_quad15g31 * v1v2;
    gpts[13][i] = _tetr_quad15g3 * v2v4 + _tetr_quad15g31 * v1v3;
    gpts[14][i] = 0.25* (v1v2 + v3v4);

  }

  weights[0]  = weights[1] = weights[ 2] = weights[ 3] = w1;
  weights[4]  = weights[5] = weights[ 6] = weights[ 7] = w2;
  weights[8]  = weights[9] = weights[10] = weights[11] = weights[12] = w3;
  weights[13] = w3;
  weights[14] = vol * _tetr_quad15w4;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the flags adapted to the given quadrature type \p qtype and the
 *         location on which the quadrature should be performed
 *
 * \param[in] qtype    \ref cs_quadrature_type_t
 * \param[in] loc      It could be \ref CS_FLAG_CELL, \ref CS_FLAG_FACE or
 *                     \ref CS_FLAG_EDGE plus \ref CS_FLAG_PRIMAL or
 *                     \ref CS_FLAG_DUAL
 *
 * \return  metadata stored in a \ref cs_eflag_t to build a \ref cs_cell_mesh_t
 */
/*----------------------------------------------------------------------------*/

cs_eflag_t
cs_quadrature_get_flag(const cs_quadrature_type_t  qtype,
                       const cs_flag_t             loc)
{
  cs_eflag_t ret_flag = 0;

  /* If necessary, enrich the mesh flag to account for the property */

  switch (qtype) {

  case CS_QUADRATURE_HIGHER:
  case CS_QUADRATURE_HIGHEST:
    ret_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ;
    /* No break, pass to the following too */
  case CS_QUADRATURE_BARY_SUBDIV:
    ret_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ;
    break;

  default:
    /* Nothing to do */
    break;

  } /* Switch */

  if (cs_flag_test(loc, cs_flag_primal_cell)) {

    switch (qtype) {

    case CS_QUADRATURE_HIGHER:
    case CS_QUADRATURE_HIGHEST:
      ret_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_HFQ;
      /* No break, pass to the following too */
    case CS_QUADRATURE_BARY_SUBDIV:
      ret_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ;
      break;

    default:
      /* Nothing to do */
      break;

    } /* Switch */

  } /* Primal cells */
  else if (cs_flag_test(loc, cs_flag_primal_face)) {

    switch (qtype) {

    case CS_QUADRATURE_HIGHER:
    case CS_QUADRATURE_HIGHEST:
      ret_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ;
      /* No break, pass to the following too */
    case CS_QUADRATURE_BARY_SUBDIV:
      ret_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ |
        CS_FLAG_COMP_PF;
      break;

    default:
      /* Nothing to do */
      break;

    } /* Switch */

  } /* Primal faces */
  else if (cs_flag_test(loc, cs_flag_dual_face)) {

    switch (qtype) {

    case CS_QUADRATURE_HIGHER:
    case CS_QUADRATURE_HIGHEST:
      ret_flag |= CS_FLAG_COMP_DFQ | CS_FLAG_COMP_SEF;
      /* No break, pass to the following too */
    case CS_QUADRATURE_BARY_SUBDIV:
      ret_flag |= CS_FLAG_COMP_EF | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_PEQ |
        CS_FLAG_COMP_PV | CS_FLAG_COMP_PE;
      break;

    default:
      /* Nothing to do */
      break;

    } /* Switch */

  } /* Dual faces */
  else if (cs_flag_test(loc, cs_flag_primal_edge)) {

    switch (qtype) {

    case CS_QUADRATURE_HIGHER:
    case CS_QUADRATURE_HIGHEST:
    case CS_QUADRATURE_BARY_SUBDIV:
      ret_flag |= CS_FLAG_COMP_PEQ;
      break;

    default:
      /* Nothing to do */
      break;

    } /* Switch */

  } /* Primal edge */

  return ret_flag;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
