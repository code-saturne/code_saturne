#ifndef __CS_CDO_SCHEME_GEOMETRY_H__
#define __CS_CDO_SCHEME_GEOMETRY_H__

/*============================================================================
 * Geometric computations for building discretization operators which is
 * shared by several files
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
 *
 * \param[in]      f_sgn    orientation of the face
 * \param[in]      pfq      quantities related to a face
 * \param[in]      deq      quantities related a dual edge
 * \param[in, out] grd_c    gradient of the Lagrange function related to xc
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_grdfc(const short int     f_sgn,
                 const cs_quant_t    pfq,
                 const cs_nvec3_t    deq,
                 cs_real_t          *grd_c)
{
  const double  hfc = cs_math_3_dot_product(pfq.unitv, deq.unitv) * deq.meas;
  const cs_real_t  ohf = -f_sgn/hfc;

  for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];
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
                  cs_real_t           *grd_v2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 *
 * \return the volume of p_{f,c}
 */
/*----------------------------------------------------------------------------*/

double
cs_compute_fwbs_q1(short int                 f,
                   const cs_cell_mesh_t     *cm,
                   cs_real_t                *wvf,
                   cs_real_t                *pefc_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] grd_c      gradient of the Lagrange function related to xc
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_fwbs_q2(short int                 f,
                   const cs_cell_mesh_t     *cm,
                   cs_real_3_t               grd_c,
                   cs_real_t                *wvf,
                   cs_real_t                *pefc_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *         wvf should be allocated to n_max_vbyc and pefc_vol to n_max_ebyf
 *
 * \param[in]      f          id of the face in the cell-wise numbering
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] grd_c      gradient of the Lagrange function related to xc
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 *
 * \return the volume of p_{f,c}
 */
/*----------------------------------------------------------------------------*/

double
cs_compute_fwbs_q3(short int                 f,
                   const cs_cell_mesh_t     *cm,
                   cs_real_3_t               grd_c,
                   cs_real_t                *wvf,
                   cs_real_t                *pefc_vol);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_SCHEME_GEOMETRY_H__ */
