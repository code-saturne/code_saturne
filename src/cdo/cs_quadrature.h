#ifndef __CS_QUADRATURE_H__
#define __CS_QUADRATURE_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_math.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_QUADRATURE_NONE,
  CS_QUADRATURE_BARY,        /* Value at the barycenter * meas */
  CS_QUADRATURE_BARY_SUBDIV, /* Value at the barycenter * meas on a sub-mesh */
  CS_QUADRATURE_HIGHER,      /* Unique weight but several Gauss points */
  CS_QUADRATURE_HIGHEST,     /* Specific weight for each Gauss points */
  CS_QUADRATURE_N_TYPES

} cs_quadrature_type_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic function pointer to compute the quadrature points for an
 *          edge from v1 -> v2
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      len      length of edge [v1, v2]
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weight (same weight for the two points)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_edge_t) (const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        double              len,
                        cs_real_3_t         gpts[],
                        double             *weights);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic functoin pointer to compute the quadrature points for a
 *          triangle
 *
 * \param[in]      v1        first vertex
 * \param[in]      v2        second vertex
 * \param[in]      v3        third vertex
 * \param[in]      area      area of triangle {v1, v2, v3}
 * \param[in, out] gpts      Gauss points
 * \param[in, out] weights   weights values
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_tria_t) (const cs_real_3_t   v1,
                        const cs_real_3_t   v2,
                        const cs_real_3_t   v3,
                        double              area,
                        cs_real_3_t         gpts[],
                        double             *weights);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to compute the quadrature points in a tetrehedra.
 *
 * \param[in]       v1       first vertex
 * \param[in]       v2       second vertex
 * \param[in]       v3       third vertex
 * \param[in]       v4       fourth vertex
 * \param[in]       vol      volume of tetrahedron {v1, v2, v3, v4}
 * \param[in, out]  gpts     Gauss points
 * \param[in, out]  weights  weigths related to each Gauss point
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_tet_t) (const cs_real_3_t  v1,
                       const cs_real_3_t  v2,
                       const cs_real_3_t  v3,
                       const cs_real_3_t  v4,
                       double             vol,
                       cs_real_3_t        gpts[],
                       double             weights[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over an edge based on a specified
 *         quadrature rule and add it to \p results
 *
 * \param[in]      tcur     current physical time of the simulation
 * \param[in]      v1       first point of the edge
 * \param[in]      v2       second point of the edge
 * \param[in]      len      length of the edge
 * \param[in]      ana      pointer to the analytic function
 * \param[in]      input    NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results  array of double
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_edge_integral_t)(double                 tcur,
                                const cs_real_3_t      v1,
                                const cs_real_3_t      v2,
                                double                 len,
                                cs_analytic_func_t    *ana,
                                void                  *input,
                                double                 results[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle based on a specified
 *         quadrature rule and add it to \p results
 *
 * \param[in]      tcur     current physical time of the simulation
 * \param[in]      v1       first point of the triangle
 * \param[in]      v2       second point of the triangle
 * \param[in]      v3       third point of the triangle
 * \param[in]      area     area of the triangle
 * \param[in]      ana      pointer to the analytic function
 * \param[in]      input    NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results  array of double
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_tria_integral_t)(double                 tcur,
                                const cs_real_3_t      v1,
                                const cs_real_3_t      v2,
                                const cs_real_3_t      v3,
                                double                 area,
                                cs_analytic_func_t    *ana,
                                void                  *input,
                                double                 results[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron based on a specified
 *         quadrature rule and add it to \p results
 *
 * \param[in]      tcur     current physical time of the simulation
 * \param[in]      v1       first point of the tetrahedron
 * \param[in]      v2       second point of the tetrahedron
 * \param[in]      v3       third point of the tetrahedron
 * \param[in]      v4       fourth point of the tetrahedron
 * \param[in]      vol      volume of the tetrahedron
 * \param[in]      ana      pointer to the analytic function
 * \param[in]      input    NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results  array of double
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_quadrature_tetra_integral_t)(double                 tcur,
                                 const cs_real_3_t      v1,
                                 const cs_real_3_t      v2,
                                 const cs_real_3_t      v3,
                                 const cs_real_3_t      v4,
                                 double                 vol,
                                 cs_analytic_func_t    *ana,
                                 void                  *input,
                                 double                 results[]);


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute constant weights for all quadratures used
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_setup(void);

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
cs_quadrature_get_type_name(const cs_quadrature_type_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for an edge from v1 -> v2 (2 points)
 *          Exact for polynomial function up to order 1
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      len      length of edge [v1, v2]
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weight (same weight for the two points)
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_edge_1pt(const cs_real_3_t  v1,
                       const cs_real_3_t  v2,
                       double             len,
                       cs_real_3_t        gpts[],
                       double            *w)
{
  gpts[0][0] = 0.5*(v1[0] + v2[0]);
  gpts[0][1] = 0.5*(v1[1] + v2[1]);
  gpts[0][2] = 0.5*(v1[2] + v2[2]);
  w[0] = len;
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
                        double            *w);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for an edge from v1 -> v2 (3 points)
 *          Exact for polynomial function up to order 5
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      len      length of edge [v1, v2]
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weights
 */
/*----------------------------------------------------------------------------*/

void
cs_quadrature_edge_3pts(const cs_real_3_t  v1,
                        const cs_real_3_t  v2,
                        double             len,
                        cs_real_3_t        gpts[],
                        double             w[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quadrature points for a triangle (1 point)
 *          Exact for polynomial function up to order 1 (barycentric approx.)
 *
 * \param[in]      v1       first vertex
 * \param[in]      v2       second vertex
 * \param[in]      v3       third vertex
 * \param[in]      area     area of triangle {v1, v2, v3}
 * \param[in, out] gpts     gauss points
 * \param[in, out] w        weight
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_1pt(const cs_real_3_t   v1,
                       const cs_real_3_t   v2,
                       const cs_real_3_t   v3,
                       double              area,
                       cs_real_3_t         gpts[],
                       double             *w)
{
  gpts[0][0] = cs_math_onethird * (v1[0] + v2[0] + v3[0]);
  gpts[0][1] = cs_math_onethird * (v1[1] + v2[1] + v3[1]);
  gpts[0][2] = cs_math_onethird * (v1[2] + v2[2] + v3[2]);
  w[0] = area;
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
                        double             *w);

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
                        double              w[]);

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
                        double              w[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the quadrature in a tetrehedra. Exact for 1st order
 *         polynomials (order 2).
 *
 * \param[in]       v1       first vertex
 * \param[in]       v2       second vertex
 * \param[in]       v3       third vertex
 * \param[in]       v4       fourth vertex
 * \param[in]       vol      volume of tetrahedron
 * \param[in, out]  gpts     1 Gauss point
 * \param[in, out]  weight   weight related to this point (= volume)
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_1pt(const cs_real_3_t   v1,
                      const cs_real_3_t   v2,
                      const cs_real_3_t   v3,
                      const cs_real_3_t   v4,
                      double              vol,
                      cs_real_3_t         gpts[],
                      double              weight[])
{
  gpts[0][0] = 0.25 * (v1[0] + v2[0] + v3[0] + v4[0]);
  gpts[0][1] = 0.25 * (v1[1] + v2[1] + v3[1] + v4[1]);
  gpts[0][2] = 0.25 * (v1[2] + v2[2] + v3[2] + v4[2]);
  weight[0] = vol;
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
                       double             weights[]);

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
                       double             weights[]);

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
                        double              weights[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over an edge with the mid-point rule and add
 *         it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the edge
 * \param[in]      v2           second point of the edge
 * \param[in]      len          length of the edge
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_edge_1pt_scal(double                 tcur,
                            const cs_real_3_t      v1,
                            const cs_real_3_t      v2,
                            double                 len,
                            cs_analytic_func_t    *ana,
                            void                  *input,
                            double                 results[])
{
  cs_real_3_t  xg;
  double  evaluation;

  // Copied from cs_quadrature_1pt
  xg[0] = .5 * (v1[0] + v2[0]);
  xg[1] = .5 * (v1[1] + v2[1]);
  xg[2] = .5 * (v1[2] + v2[2]);

  ana(tcur, 1, NULL, xg, false, input, &evaluation);

  *results += len * evaluation;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over an edge with a quadrature rule using 2
 *         Gauss points and a unique weight and add it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the edge
 * \param[in]      v2           second point of the edge
 * \param[in]      len          length of the edge
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_edge_2pts_scal(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             double                len,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[2];
  double  evaluation[2], weights[2];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_edge_2pts(v1, v2, len, gauss_pts, weights);

  ana(tcur, 2, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  /* Return results */
  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over an edge with a quadrature rule using 3
 *         Gauss points and weights and add it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the edge
 * \param[in]      v2           second point of the edge
 * \param[in]      len          length of the edge
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_edge_3pts_scal(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             double                len,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[3];
  double  evaluation[3], weights[3];

  /* Compute Gauss points and its weights */
  cs_quadrature_edge_3pts(v1, v2, len, gauss_pts, weights);

  ana(tcur, 3, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1] +
              weights[2] * evaluation[2];
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle using a barycentric quadrature
 *         rule and add it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_1pt_scal(double                 tcur,
                            const cs_real_3_t      v1,
                            const cs_real_3_t      v2,
                            const cs_real_3_t      v3,
                            double                 area,
                            cs_analytic_func_t    *ana,
                            void                  *input,
                            double                 results[])
{
  cs_real_3_t  xg;
  double  evaluation;

  /* Copied from cs_quadrature_1pt */
  xg[0] = cs_math_onethird * (v1[0] + v2[0] + v3[0]);
  xg[1] = cs_math_onethird * (v1[1] + v2[1] + v3[1]);
  xg[2] = cs_math_onethird * (v1[2] + v2[2] + v3[2]);

  ana(tcur, 1, NULL, xg, false, input, &evaluation);

  *results += area * evaluation;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         3 Gauss points and a unique weight and add it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_3pts_scal(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[3];
  double  evaluation[3], weights[3];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tria_3pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 3, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  /* Return results */
  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1] +
              weights[2] * evaluation[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         4 Gauss points and 4 weights and add it to \p results
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_4pts_scal(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[4], weights[4];

  /* Compute Gauss points and its weights */
  cs_quadrature_tria_4pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1] +
              weights[2] * evaluation[2] + weights[3] * evaluation[3];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle using a barycentric quadrature
 *         rule and add it to \p results
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results  array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_1pt_vect(double                 tcur,
                            const cs_real_3_t      v1,
                            const cs_real_3_t      v2,
                            const cs_real_3_t      v3,
                            double                 area,
                            cs_analytic_func_t    *ana,
                            void                  *input,
                            double                 results[])
{
  cs_real_3_t  xg;
  double evaluation[3];

  /* Copied from cs_quadrature_1pt */
  xg[0] = cs_math_onethird * (v1[0] + v2[0] + v3[0]);
  xg[1] = cs_math_onethird * (v1[1] + v2[1] + v3[1]);
  xg[2] = cs_math_onethird * (v1[2] + v2[2] + v3[2]);

  ana(tcur, 1, NULL, xg, false, input, evaluation);

  results[0] += area * evaluation[0];
  results[1] += area * evaluation[1];
  results[2] += area * evaluation[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         3 Gauss points and a unique weight and add it to \p results
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_3pts_vect(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[3];
  double  evaluation[3*3], weights[3];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tria_3pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 3, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 3; p++) {
    results[0] += weights[p] * evaluation[3*p  ];
    results[1] += weights[p] * evaluation[3*p+1];
    results[2] += weights[p] * evaluation[3*p+2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         4 Gauss points and 4 weights and add it to \p results.
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_4pts_vect(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[3*4], weights[4];

  /* Compute Gauss points and its weights */
  cs_quadrature_tria_4pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 4; p++) {
    results[0] += weights[p] * evaluation[3*p  ];
    results[1] += weights[p] * evaluation[3*p+1];
    results[2] += weights[p] * evaluation[3*p+2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle using a barycentric quadrature
 *         rule and add it to \p results.
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_1pt_tens(double                 tcur,
                            const cs_real_3_t      v1,
                            const cs_real_3_t      v2,
                            const cs_real_3_t      v3,
                            double                 area,
                            cs_analytic_func_t    *ana,
                            void                  *input,
                            double                 results[])
{
  cs_real_3_t  xg;
  double evaluation[9];

  /* Copied from cs_quadrature_1pt */
  xg[0] = cs_math_onethird * (v1[0] + v2[0] + v3[0]);
  xg[1] = cs_math_onethird * (v1[1] + v2[1] + v3[1]);
  xg[2] = cs_math_onethird * (v1[2] + v2[2] + v3[2]);

  ana(tcur, 1, NULL, xg, false, input, evaluation);

  for (short int ij = 0; ij < 9; ij++)
    results[ij] += area * evaluation[ij];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         3 Gauss points and a unique weight and add it to \p results.
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_3pts_tens(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[3];
  double  evaluation[9*3], weights[3];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tria_3pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 3, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 3; p++) {
    const double wp = weights[p];
    double *eval_p = evaluation + 9*p;
    for (short int ij = 0; ij < 9; ij++)
      results[ij] += wp * eval_p[ij];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a triangle with a quadrature rule using
 *         4 Gauss points and 4 weights and add it to \p results.
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the triangle
 * \param[in]      v2           second point of the triangle
 * \param[in]      v3           third point of the triangle
 * \param[in]      area         area of the triangle
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tria_4pts_tens(double                tcur,
                             const cs_real_3_t     v1,
                             const cs_real_3_t     v2,
                             const cs_real_3_t     v3,
                             double                area,
                             cs_analytic_func_t   *ana,
                             void                 *input,
                             double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[9*4], weights[4];

  /* Compute Gauss points and its weights */
  cs_quadrature_tria_4pts(v1, v2, v3, area, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 4; p++) {
    const double wp = weights[p];
    double *eval_p = evaluation + 9*p;
    for (short int ij = 0; ij < 9; ij++)
      results[ij] += wp * eval_p[ij];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron using a barycentric
 *         quadrature rule and add it to \p results.
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_1pt_scal(double                 tcur,
                           const cs_real_3_t      v1,
                           const cs_real_3_t      v2,
                           const cs_real_3_t      v3,
                           const cs_real_3_t      v4,
                           double                 vol,
                           cs_analytic_func_t    *ana,
                           void                  *input,
                           double                 results[])
{
  cs_real_3_t  xg;
  double  evaluation;

  /* Copied from cs_quadrature_tet_1pt */
  xg[0] = 0.25 * (v1[0] + v2[0] + v3[0] + v4[0]);
  xg[1] = 0.25 * (v1[1] + v2[1] + v3[1] + v4[1]);
  xg[2] = 0.25 * (v1[2] + v2[2] + v3[2] + v4[2]);

  ana(tcur, 1, NULL, xg, false, input, &evaluation);

  *results += vol * evaluation;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         4 Gauss points and a unique weight and add it to \p results.
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_4pts_scal(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[4], weights[4];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_4pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  /* Return results */
  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1] +
              weights[2] * evaluation[2] + weights[3] * evaluation[3];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         5 Gauss points and 5 weights and add it to \p results.
 *         Case of a scalar-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_5pts_scal(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[5];
  double  evaluation[5], weights[5];

  /* Compute Gauss points and its weights */
  cs_quadrature_tet_5pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 5, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  *results += weights[0] * evaluation[0] + weights[1] * evaluation[1] +
              weights[2] * evaluation[2] + weights[3] * evaluation[3] +
              weights[4] * evaluation[4];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron using a barycentric
 *         quadrature rule and add it to \p results.
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_1pt_vect(double                 tcur,
                           const cs_real_3_t      v1,
                           const cs_real_3_t      v2,
                           const cs_real_3_t      v3,
                           const cs_real_3_t      v4,
                           double                 vol,
                           cs_analytic_func_t    *ana,
                           void                  *input,
                           double                 results[])
{
  cs_real_3_t  xg;
  double evaluation[3];

  /* Copied from cs_quadrature_tet_1pt */
  xg[0] = 0.25 * (v1[0] + v2[0] + v3[0] + v4[0]);
  xg[1] = 0.25 * (v1[1] + v2[1] + v3[1] + v4[1]);
  xg[2] = 0.25 * (v1[2] + v2[2] + v3[2] + v4[2]);

  ana(tcur, 1, NULL, xg, false, input, evaluation);

  results[0] += vol * evaluation[0];
  results[1] += vol * evaluation[1];
  results[2] += vol * evaluation[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         4 Gauss points and a unique weight and add it to \p results.
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_4pts_vect(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[3*4], weights[4];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_4pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 4; p++) {
    results[0] += weights[p] * evaluation[3*p  ];
    results[1] += weights[p] * evaluation[3*p+1];
    results[2] += weights[p] * evaluation[3*p+2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         5 Gauss points and 5 weights and add it to \p results.
 *         Case of a vector-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_5pts_vect(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[5];
  double  evaluation[3*5], weights[5];

  /* Compute Gauss points and its weights */
  cs_quadrature_tet_5pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 5, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 5; p++) {
    results[0] += weights[p] * evaluation[3*p  ];
    results[1] += weights[p] * evaluation[3*p+1];
    results[2] += weights[p] * evaluation[3*p+2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron using a barycentric
 *         quadrature rule and add it to \p results.
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_1pt_tens(double                 tcur,
                           const cs_real_3_t      v1,
                           const cs_real_3_t      v2,
                           const cs_real_3_t      v3,
                           const cs_real_3_t      v4,
                           double                 vol,
                           cs_analytic_func_t    *ana,
                           void                  *input,
                           double                 results[])
{
  cs_real_3_t  xg;
  double evaluation[9];

  /* Copied from cs_quadrature_tet_1pt */
  xg[0] = 0.25 * (v1[0] + v2[0] + v3[0] + v4[0]);
  xg[1] = 0.25 * (v1[1] + v2[1] + v3[1] + v4[1]);
  xg[2] = 0.25 * (v1[2] + v2[2] + v3[2] + v4[2]);

  ana(tcur, 1, NULL, xg, false, input, evaluation);

  for (short int ij = 0; ij < 9; ij++)
    results[ij] += vol * evaluation[ij];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         4 Gauss points and a unique weight and add it to \p results.
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_4pts_tens(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[4];
  double  evaluation[9*4], weights[4];

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_4pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 4, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 4; p++) {
    const double wp = weights[p];
    double *eval_p = evaluation + 9*p;
    for (short int ij = 0; ij < 9; ij++)
      results[ij] += wp * eval_p[ij];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron with a quadrature rule using
 *         5 Gauss points and 5 weights and add it to \p results
 *         Case of a tensor-valued function.
 *
 * \param[in]      tcur         current physical time of the simulation
 * \param[in]      v1           first point of the tetrahedron
 * \param[in]      v2           second point of the tetrahedron
 * \param[in]      v3           third point of the tetrahedron
 * \param[in]      v4           fourth point of the tetrahedron
 * \param[in]      vol          volume of the tetrahedron
 * \param[in]      ana          pointer to the analytic function
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] results      array of double
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_quadrature_tet_5pts_tens(double                tcur,
                            const cs_real_3_t     v1,
                            const cs_real_3_t     v2,
                            const cs_real_3_t     v3,
                            const cs_real_3_t     v4,
                            double                vol,
                            cs_analytic_func_t   *ana,
                            void                 *input,
                            double                results[])
{
  cs_real_3_t  gauss_pts[5];
  double  evaluation[9*5], weights[5];

  /* Compute Gauss points and its weights */
  cs_quadrature_tet_5pts(v1, v2, v3, v4, vol, gauss_pts, weights);

  ana(tcur, 5, NULL, (const cs_real_t *)gauss_pts, false, input, evaluation);

  for (int p = 0; p < 5; p++) {
    const double wp = weights[p];
    double *eval_p = evaluation + 9*p;
    for (short int ij = 0; ij < 9; ij++)
      results[ij] += wp * eval_p[ij];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the integral function according to the the quadrature type
 *         and the stride provided
 *
 * \param[in]     dim          dimension of the function to integrate
 * \param[in]     qtype        quadrature type
 *
 * \return a pointer to the integral function
 */
/*----------------------------------------------------------------------------*/

static inline cs_quadrature_tria_integral_t *
cs_quadrature_get_integral_function_tria(short int             dim,
                                         cs_quadrature_type_t  qtype)
{
  switch (dim) {

  case 1: /* Scalar-valued integral */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tria_1pt_scal;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tria_3pts_scal;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tria_4pts_scal;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  case 3: /* Vector-valued case */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tria_1pt_vect;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tria_3pts_vect;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tria_4pts_vect;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  case 9: /* Tensor-valued case */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tria_1pt_tens;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tria_3pts_tens;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tria_4pts_tens;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid dimension value %d. Only 1, 3 and 9 are valid.\n",
              __func__, dim);

  } /* switch on dim */

  return NULL; /* Should not go to this stage */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the integral function according to the the quadrature type
 *         and the stride provided
 *
 * \param[in]     dim          dimension of the function to integrate
 * \param[in]     qtype        quadrature type
 *
 * \return a pointer to the integral function
 */
/*----------------------------------------------------------------------------*/

static inline cs_quadrature_tetra_integral_t*
cs_quadrature_get_integral_function_tetra(short int             dim,
                                          cs_quadrature_type_t  qtype)
{
  switch (dim) {

  case 1: /* Scalar-valued case */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tet_1pt_scal;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tet_4pts_scal;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tet_5pts_scal;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  case 3: /* Vector-valued case */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tet_1pt_vect;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tet_4pts_vect;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tet_5pts_vect;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  case 9: /* Tensor-valued case */

    switch (qtype) {

    case CS_QUADRATURE_BARY:
    case CS_QUADRATURE_BARY_SUBDIV:
      return cs_quadrature_tet_1pt_tens;
    case CS_QUADRATURE_HIGHER:
      return cs_quadrature_tet_4pts_tens;
    case CS_QUADRATURE_HIGHEST:
      return cs_quadrature_tet_5pts_tens;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid quadrature type\n", __func__);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid dimension value %d. Only 1, 3 and 9 are valid.\n",
              __func__, dim);

  } /* Switch on dim */

  /* Avoid no return warning */
  return NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_QUADRATURE_H__ */
