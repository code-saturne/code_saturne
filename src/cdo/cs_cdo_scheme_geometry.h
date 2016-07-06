#ifndef __CS_CDO_SCHEME_GEOMETRY_H__
#define __CS_CDO_SCHEME_GEOMETRY_H__

/*============================================================================
 * Geometric computations for building discretization operators which is
 * shared by several files
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the constant gradient of the Lagrange function
 *         attached to xc in p_{f,c} (constant inside this volume)
 *
 * \param[in]      fm       pointer to a cs_face_mesh_t structure
 * \param[in, out] grd_c    gradient of the Lagrange function related to xc
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_compute_grdc(const cs_face_mesh_t     *fm,
                cs_real_t                *grd_c)
{
  const cs_quant_t  pfq = fm->face;
  const cs_nvec3_t  deq = fm->dedge;
  const double  hf = cs_math_3_dot_product(pfq.unitv, deq.unitv) * deq.meas;
  const cs_real_t  ohf = -fm->f_sgn/hf;

  for (int k = 0; k < 3; k++)
    grd_c[k] = ohf * pfq.unitv[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute tef (the area of the triangle of base e and apex f)
 *
 * \param[in]  e     id of the edge in the face-wise numbering
 * \param[in]  fm    pointer to a cs_face_mesh_t structure
 *
 * \return the value the area of the triangle tef
 */
/*----------------------------------------------------------------------------*/

inline static double
cs_compute_tef(short int                 e,
               const cs_face_mesh_t     *fm)
{
  double  xef_len;
  cs_real_3_t  xef_un, un;

  const cs_quant_t  pfq = fm->face;
  const cs_quant_t  peq = fm->edge[e];

  cs_math_3_length_unitv(peq.center, pfq.center, &xef_len, xef_un);
  cs_math_3_cross_product(xef_un, peq.unitv, un);

  /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
  return 0.5 * xef_len * peq.meas * cs_math_3_norm(un);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute tec (the area of the triangle of base e and apex c)
 *
 * \param[in]  e     id of the edge in the face-wise numbering
 * \param[in]  cm    pointer to a cs_cell_mesh_t structure
 *
 * \return the value the area of the triangle tec
 */
/*----------------------------------------------------------------------------*/

inline static double
cs_compute_tec(short int                 e,
               const cs_cell_mesh_t     *cm)
{
  double  xec_len;
  cs_real_3_t  xec_un, un;

  const cs_quant_t  peq = cm->edge[e];

  cs_math_3_length_unitv(peq.center, cm->xc, &xec_len, xec_un);
  cs_math_3_cross_product(xec_un, peq.unitv, un);

  /* tec = ||(xe -xc) x e||/2 */
  return 0.5 * xec_len * peq.meas * cs_math_3_norm(un);
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
                  cs_real_t           *grd_v2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute tef (the are of the triangle of base e and apex f
 *         Compute also the value of the constant gradient attached to xc in
 *         p_{f,c}
 *
 * \param[in]      f        id of the face in the cell-wise numbering
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out] tef      pointer to an array storing area of tef triangles
 * \param[in, out] grd_c    gradient of the Lagrange function related to xc
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_tef_grdc(short int                 f,
                    const cs_cell_mesh_t     *cm,
                    cs_real_t                *tef,
                    cs_real_t                *grd_c);

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
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_fwbs_q0(short int                 f,
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
