#ifndef __CS_CDO_SCHEME_GEOMETRY_H__
#define __CS_CDO_SCHEME_GEOMETRY_H__

/*============================================================================
 * Geometric computations for building discretization operators which is
 * shared by several files
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Inline static function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the constant gradient of the Lagrange function
 *         attached to xc in p_{f,c} (constant inside this volume)
 *         Cellwise version.
 *
 * \param[in]      f        face number in the cellwise numbering to handle
 * \param[in]      cm       pointer to a cell_mesh_t structure
 * \param[in, out] grd_c    gradient of the Lagrange function related to xc
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_grdfc_cw(short int               f,
                    const cs_cell_mesh_t   *cm,
                    cs_real_t              *grd_c)
{
  /* Sanity checks */
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_HFQ));

  const cs_real_t  ohf = -cm->f_sgn[f]/cm->hfc[f];

  grd_c[0] = ohf * cm->face[f].unitv[0];
  grd_c[1] = ohf * cm->face[f].unitv[1];
  grd_c[2] = ohf * cm->face[f].unitv[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the constant gradient of the Lagrange function
 *         attached to xc in p_{f,c} (constant inside this volume)
 *         Facewise version.
 *
 * \param[in]      f        face number in the cellwise numbering to handle
 * \param[in]      cm       pointer to a cell_mesh_t structure
 * \param[in, out] grd_c    gradient of the Lagrange function related to xc
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_grdfc_fw(const cs_face_mesh_t   *fm,
                    cs_real_t              *grd_c)
{
  const cs_real_t  ohf = -fm->f_sgn/fm->hfc;

  grd_c[0] = ohf * fm->face.unitv[0];
  grd_c[1] = ohf * fm->face.unitv[1];
  grd_c[2] = ohf * fm->face.unitv[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weight wef = |tef|/|f|
 *
 * \param[in]       f          id of the face in the cell-wise numbering
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out]  wef        pointer to an array storing the weight/vertex
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_wef(short int                 f,
               const cs_cell_mesh_t     *cm,
               cs_real_t                *wef)
{
  /* Sanity checks */
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FEQ));

  const short int  *f2e_idx = cm->f2e_idx + f;
  const double  *tef_vals = cm->tef + f2e_idx[0];
  const double  inv_f = 1./cm->face[f].meas;

  /* Compute a weight for each edge of the current face */
  for (short int e = 0; e < f2e_idx[1] - f2e_idx[0]; e++)
    wef[e] = tef_vals[e] * inv_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume related to each tetrahedron of base tef and apex
 *         x_c (there are as many tetrahedra as edges in a face)
 *
 * \param[in]       f      id of the face in the cell-wise numbering
 * \param[in]       cm     pointer to a cs_cell_mesh_t structure
 * \param[in, out]  pefc   pointer to an array storing the volumes
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_pefc(short int                 f,
                const cs_cell_mesh_t     *cm,
                cs_real_t                *pefc)
{
  /* Sanity checks */
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_PFC));

  const short int  *f2e_idx = cm->f2e_idx + f;
  const double  *tef_vals = cm->tef + f2e_idx[0];
  const double  f_coef = cm->pvol_f[f]/cm->face[f].meas;

  /* Compute a weight for each edge of the current face */
  for (short int e = 0; e < f2e_idx[1] - f2e_idx[0]; e++)
    pefc[e] = tef_vals[e] * f_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *
 * \param[in]       f          id of the face in the cell-wise numbering
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out]  wvf        pointer to an array storing the weight/vertex
 *                             (allocated to the number of cell vertices)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_wvf(short int                 f,
               const cs_cell_mesh_t     *cm,
               cs_real_t                *wvf)
{
  /* Sanity checks */
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  /* Reset weights */
  memset(wvf, 0, cm->n_vc*sizeof(cs_real_t));

  const short int  *f2e_idx = cm->f2e_idx + f;
  const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
  const double  *tef_vals = cm->tef + f2e_idx[0];
  const double  inv_f = 1./cm->face[f].meas;

  /* Compute a weight for each vertex of the current face */
  for (short int e = 0; e < f2e_idx[1] - f2e_idx[0]; e++) {

    const short int  *v = cm->e2v_ids + 2*f2e_ids[e];
    const double  w = 0.5*tef_vals[e] * inv_f;
    wvf[v[0]] += w;
    wvf[v[1]] += w;

  }
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
                                  cs_real_t               cov[3]);

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
                          cs_real_t               inertia[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the gradient of a Lagrange function related to primal
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
                  cs_real_t           *grd_v2);

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
                   cs_real_t                *wef);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_SCHEME_GEOMETRY_H__ */
