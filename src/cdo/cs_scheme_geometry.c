/*============================================================================
 * Geometric computations for building discretization operators which is
 * shared by several files
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_quadrature.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_scheme_geometry.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_SCHEME_GEOMETRY_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the covariance tensor with the contribution of the current
 *         triangle
 *
 * \param[in]      x1       1st vertex coordinate
 * \param[in]      x2       2nd vertex coordinate
 * \param[in]      x3       3rd vertex coordinate
 * \param[in]      ax       main X-axis for the face-related coordinate system
 * \param[in]      ay       main Y-axis for the face-related coordinate system
 * \param[in]      center   center used for the computation
 * \param[in]      area     area of the triangle
 * \param[in, out] tensor   covariance tensor to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tria_to_covariance(const cs_real_t     x1[3],
                        const cs_real_t     x2[3],
                        const cs_real_t     x3[3],
                        const cs_nvec3_t    ax,
                        const cs_nvec3_t    ay,
                        const cs_real_t     center[3],
                        cs_real_t           area,
                        cs_real_t           tensor[3])
{
  cs_real_3_t gpts[3], r;
  cs_real_t   gw;

  cs_quadrature_tria_3pts(x1, x2, x3, area,  /* 2nd-order exact */
                          gpts, &gw);

  for (short int gp = 0; gp < 3; gp++) {

    for (int k = 0; k < 3; k++) r[k] = gpts[gp][k] - center[k];

    const cs_real_t  xf = ax.meas * _dp3(ax.unitv, r);
    const cs_real_t  yf = ay.meas * _dp3(ay.unitv, r);

    tensor[0] += gw * xf*xf;
    tensor[1] += gw * xf*yf;
    tensor[2] += gw * yf*yf;

  } /* Loop on gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the computation of the inertia tensor with the contribution
 *         of a tetrahedron
 *
 * \param[in]      x1       1st vertex coordinate
 * \param[in]      x2       2nd vertex coordinate
 * \param[in]      x3       3rd vertex coordinate
 * \param[in]      x4       4th vertex coordinate
 * \param[in]      center   center used for the computation
 * \param[in]      vol      volume of the tetrahedron
 * \param[in, out] tensor   inertia tensor to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tetra_to_inertia2(const cs_real_t    x1[3],
                       const cs_real_t    x2[3],
                       const cs_real_t    x3[3],
                       const cs_real_t    x4[3],
                       const cs_real_t    center[3],
                       cs_real_t          vol,
                       cs_real_33_t       tensor)
{
  cs_real_3_t gpts[4], r;
  cs_real_t   _gw[4];

  cs_quadrature_tet_4pts(x1, x2, x3, x4, vol, gpts, _gw);

  const cs_real_t gw = _gw[0];  /* same weight for all Gauss points */
  for (short int gp = 0; gp < 4; gp++) {

    for (int k = 0; k < 3; k++) r[k] = gpts[gp][k] - center[k];

    const cs_real_t  rx2 = r[0]*r[0], ry2 = r[1]*r[1], rz2 = r[2]*r[2];

    tensor[0][0] += gw * (ry2 + rz2), tensor[0][1] -= gw * r[0]*r[1];
    tensor[0][2] -= gw * r[0]*r[2],   tensor[1][1] += gw * (rx2 + rz2);
    tensor[1][2] -= gw * r[1]*r[2],   tensor[2][2] += gw * (rx2 + ry2);

  } /* Loop on gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the computation of the inertia tensor with the contribution
 *         of a tetrahedron
 *
 * \param[in]      x1       1st vertex coordinate
 * \param[in]      x2       2nd vertex coordinate
 * \param[in]      x3       3rd vertex coordinate
 * \param[in]      x4       4th vertex coordinate
 * \param[in]      center   center used for the computation
 * \param[in]      vol      volume of the tetrahedron
 * \param[in, out] tensor   inertia tensor to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tetra_to_inertia3(const cs_real_t    x1[3],
                       const cs_real_t    x2[3],
                       const cs_real_t    x3[3],
                       const cs_real_t    x4[3],
                       const cs_real_t    center[3],
                       cs_real_t          vol,
                       cs_real_33_t       tensor)
{
  cs_real_3_t gpts[4], r;
  cs_real_t   _gw[4];

  cs_quadrature_tet_4pts(x1, x2, x3, x4, vol, gpts, _gw);

  const cs_real_t gw = _gw[0];  /* same weight for all Gauss points */
  for (short int gp = 0; gp < 4; gp++) {

    for (int k = 0; k < 3; k++) r[k] = gpts[gp][k] - center[k];

    tensor[0][0] += gw * r[0]*r[0], tensor[0][1] += gw * r[0]*r[1];
    tensor[0][2] += gw * r[0]*r[2], tensor[1][1] += gw * r[1]*r[1];
    tensor[1][2] += gw * r[1]*r[2], tensor[2][2] += gw * r[2]*r[2];

  } /* Loop on gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the computation of the inertia tensor with the contribution
 *         of a tetrahedron
 *         Ref.: F. Tonon "Explicit exact formulas for the 3D tetrahedron
 *         inertia tensor in terms of its vertex coordinates" (2004)
 *         J. of Mathematics and Statistics
 *
 * \param[in]      x1       1st vertex coordinate
 * \param[in]      x2       2nd vertex coordinate
 * \param[in]      x3       3rd vertex coordinate
 * \param[in]      x4       4th vertex coordinate
 * \param[in]      center   center used for the computation
 * \param[in]      vol      volume of the tetrahedron
 * \param[in, out] tensor   inertia tensor to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tetra_to_inertia(const cs_real_t    x1[3],
                      const cs_real_t    x2[3],
                      const cs_real_t    x3[3],
                      const cs_real_t    x4[3],
                      const cs_real_t    center[3],
                      cs_real_t          vol,
                      cs_real_33_t       tensor)
{
  const cs_real_t  vol_weight = 0.1 * vol, vol_wby2 = 0.5 * vol_weight;

  const cs_real_3_t  r1 = {x1[0]-center[0], x1[1]-center[1], x1[2]-center[2]};
  const cs_real_3_t  r2 = {x2[0]-center[0], x2[1]-center[1], x2[2]-center[2]};
  const cs_real_3_t  r3 = {x3[0]-center[0], x3[1]-center[1], x3[2]-center[2]};
  const cs_real_3_t  r4 = {x4[0]-center[0], x4[1]-center[1], x4[2]-center[2]};

  const cs_real_t  x_2 = r1[0]*r1[0] + r2[0]*r2[0] + r3[0]*r3[0] + r4[0]*r4[0];
  const cs_real_t  y_2 = r1[1]*r1[1] + r2[1]*r2[1] + r3[1]*r3[1] + r4[1]*r4[1];
  const cs_real_t  z_2 = r1[2]*r1[2] + r2[2]*r2[2] + r3[2]*r3[2] + r4[2]*r4[2];

  const cs_real_t  xx =
    r1[0]*(r2[0] + r3[0] + r4[0]) + r2[0]*(r3[0] + r4[0]) + r3[0]*r4[0];
  const cs_real_t  yy =
    r1[1]*(r2[1] + r3[1] + r4[1]) + r2[1]*(r3[1] + r4[1]) + r3[1]*r4[1];
  const cs_real_t  zz =
    r1[2]*(r2[2] + r3[2] + r4[2]) + r2[2]*(r3[2] + r4[2]) + r3[2]*r4[2];

  tensor[0][0] += vol_weight * (y_2 + yy + z_2 + zz);
  tensor[1][1] += vol_weight * (x_2 + xx + z_2 + zz);
  tensor[2][2] += vol_weight * (x_2 + xx + y_2 + yy);

  tensor[1][2] += vol_wby2 * ( r1[1] * (2*r1[2] +   r2[2] +   r3[2] +   r4[2]) +
                               r2[1] * (  r1[2] + 2*r2[2] +   r3[2] +   r4[2]) +
                               r3[1] * (  r1[2] +   r2[2] + 2*r3[2] +   r4[2]) +
                               r4[1] * (  r1[2] +   r2[2] +   r3[2] + 2*r4[2]) );
  tensor[0][2] += vol_wby2 * ( r1[0] * (2*r1[2] +   r2[2] +   r3[2] +   r4[2]) +
                               r2[0] * (  r1[2] + 2*r2[2] +   r3[2] +   r4[2]) +
                               r3[0] * (  r1[2] +   r2[2] + 2*r3[2] +   r4[2]) +
                               r4[0] * (  r1[2] +   r2[2] +   r3[2] + 2*r4[2]) );
  tensor[0][1] += vol_wby2 * ( r1[0] * (2*r1[1] +   r2[1] +   r3[1] +   r4[1]) +
                               r2[0] * (  r1[1] + 2*r2[1] +   r3[1] +   r4[1]) +
                               r3[0] * (  r1[1] +   r2[1] + 2*r3[1] +   r4[1]) +
                               r4[0] * (  r1[1] +   r2[1] +   r3[1] + 2*r4[1]) );
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the inertial matrix of a cell with respect to the point
 *          called "center". This computation is performed exactly thanks to
 *          quadrature based on a "tetrahedrization" of the cell.
 *
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       f        id of the face in the cell numbering
 * \param[in]       ax       main X-axis for the face-related coordinate system
 * \param[in]       ay       main Y-axis for the face-related coordinate system
 * \param[in]       center   coordinates of the face center
 * \param[in, out]  cov      2x2 symmetric covariance matrix to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_face_covariance_tensor(const cs_cell_mesh_t   *cm,
                                  short int               f,
                                  const cs_nvec3_t        ax,
                                  const cs_nvec3_t        ay,
                                  const cs_real_t         center[3],
                                  cs_real_t               cov[3])
{
  assert(cm != NULL && f > -1); /* Sanity checks */

  cov[0] = cov[1] = cov[2] = 0;

  const int  start = cm->f2e_idx[f];
  const int  end = cm->f2e_idx[f+1];
  const short int  n_ef = end - start; /* ertices (=#edges) */
  const short int  *f2e_ids = cm->f2e_ids + start;
  const cs_quant_t  pfq = cm->face[f];

  /* Switching on face-type: optimized version for triangles */
  switch (n_ef) {
  case CS_TRIANGLE_CASE: /* Triangle */
    {
      short int  v0, v1, v2;
      cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

      _add_tria_to_covariance(cm->xv + 3*v0, cm->xv +3*v1, cm->xv + 3*v2,
                              ax, ay, center, pfq.meas,
                              cov);
    }
    break;

  default:
    {
      assert(n_ef > 3);
      const double  *tef = cm->tef + start;

      for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

        /* Edge-related variables */
        const short int e0 = f2e_ids[e];
        const short int v0 = cm->e2v_ids[2*e0];
        const short int v1 = cm->e2v_ids[2*e0+1];

        _add_tria_to_covariance(cm->xv + 3*v0, cm->xv +3*v1, pfq.center,
                                ax, ay, center, tef[e],
                                cov);

      }
    }
    break;

  } /* End of switch on n_ef */

#if defined(DEBUG) && !defined(NDEBUG) && CS_SCHEME_GEOMETRY_DBG > 0
  printf("\n           |% 5.3e % 5.3e|\n"
         "   cov(%2d) |% 5.3e % 5.3e|\n",
         cov[0], cov[1], f, cov[1], cov[2]);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the inertial matrix of a cell with respect to the point
 *          called "center". This computation is performed exactly thanks to
 *          quadrature based on a "tetrahedrization" of the cell.
 *
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       center   coordinates of the cell center
 * \param[in, out]  inertia  inertia matrix to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_inertia_tensor(const cs_cell_mesh_t   *cm,
                          const cs_real_t         center[3],
                          cs_real_t               inertia[3][3])
{
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PV | CS_FLAG_COMP_PFQ |
                       CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ));

  cs_real_33_t  M = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  /* Switching on cell-type: optimised version for tetra */

  switch (cm->type) {

  case FVM_CELL_TETRA:
    _add_tetra_to_inertia3(cm->xv, cm->xv+3, cm->xv+6, cm->xv+9,
                           center, cm->vol_c,
                           M);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_1ov3 * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; /* #vertices (=#edges) */
      const short int *f2e_ids = cm->f2e_ids + start;
      assert(n_vf > 2);

      switch(n_vf){

      case CS_TRIANGLE_CASE: /* Triangle */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _add_tetra_to_inertia3(cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, cm->xc,
                                 center, hf_coef * pfq.meas,
                                 M);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            /* Edge-related variables */
            const short int e0  = f2e_ids[e];
            const short int v0 = cm->e2v_ids[2*e0];
            const short int v1 = cm->e2v_ids[2*e0+1];

            _add_tetra_to_inertia3(cm->xv+3*v0, cm->xv+3*v1, pfq.center, cm->xc,
                                   center, hf_coef * tef[e],
                                   M);

          }
        }
        break;

      } /* End of switch */
    }   /* End of loop on faces */

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */

  /* Inertia tensor is symmetric */
  for (short int i = 0; i < 3; ++i) {
    inertia[i][i] = M[i][i];
    for (short int j = 0; j < i; ++j)
      inertia[i][j] = inertia[j][i] = M[j][i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the gradient of a Lagrange hat function related to primal
 *          vertices in a p_{ef,c} subvolume of a cell c where e is an edge
 *          belonging to the face f with vertices v1 and v2
 *
 * \param[in]       v1        number of the first vertex in cell numbering
 * \param[in]       v2        number of the second vertex in cell numbering
 * \param[in]       deq       dual edge quantities
 * \param[in]       uvc       xc --> xv unit tangent vector
 * \param[in]       lvc       xc --> xv vector length
 * \param[in, out]  grd_v1   gradient of Lagrange function related to v1
 * \param[in, out]  grd_v2   gradient of Lagrange function related to v2
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_grd_ve(const short int      v1,
                  const short int      v2,
                  const cs_nvec3_t     deq,
                  const cs_real_3_t    uvc[],
                  const cs_real_t      lvc[],
                  cs_real_t           *grd_v1,
                  cs_real_t           *grd_v2)
{
  double  hv;
  cs_real_3_t  unormal;

  /* Gradient for v1
     Normal direction to the plane in opposition to v1
     Height from this plane to the vertex v1 */
  cs_math_3_cross_product(uvc[v2], deq.unitv, unormal);
  hv = lvc[v1] * _dp3(uvc[v1], unormal);
  assert(fabs(hv) > DBL_EPSILON); /* Sanity check */

  const double  ohv1 = 1/hv;
  for (int k = 0; k < 3; k++) grd_v1[k] = unormal[k] * ohv1;

  /* Gradient for v2
     Normal direction to the plane in opposition to v2
     Height from this plane to the vertex v2 */
  cs_math_3_cross_product(uvc[v1], deq.unitv, unormal);
  hv = lvc[v2] * _dp3(uvc[v2], unormal);
  assert(fabs(hv) > DBL_EPSILON); /* Sanity check */

  const double  ohv2 = 1/hv;
  for (int k = 0; k < 3; k++) grd_v2[k] = unormal[k] * ohv2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f} and
 *         the weight related to each edge
 *         w_{v,f} = |dc(v) cap f|/|f|
 *         Sum of w_{v,f} over the face vertices is equal to 1
 *         Sum of w_{e,f} over the face edges is equal to 1
 *
 * \param[in]       f      id of the face in the cell-wise numbering
 * \param[in]       cm     pointer to a cs_cell_mesh_t structure
 * \param[in, out]  wvf    weights of each face vertex
 * \param[in, out]  wef    weights of each face edge
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_wef_wvf(short int                 f,
                   const cs_cell_mesh_t     *cm,
                   cs_real_t                *wvf,
                   cs_real_t                *wef)
{
  /* Reset weights */
  memset(wvf, 0, cm->n_vc*sizeof(cs_real_t));
  memset(wef, 0, cm->n_ec*sizeof(cs_real_t));

  const short int  *f2e_idx = cm->f2e_idx + f;
  const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
  const double  *tef_vals = cm->tef + f2e_idx[0];
  const double  inv_f = 1./cm->face[f].meas;

  /* Compute a weight for each vertex of the current face */
  for (short int e = 0; e < f2e_idx[1] - f2e_idx[0]; e++) {

    const short int  *v = cm->e2v_ids + 2*f2e_ids[e];
    wef[e] = tef_vals[e] * inv_f;
    wvf[v[0]] += 0.5*wef[e];
    wvf[v[1]] += 0.5*wef[e];

  }
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
