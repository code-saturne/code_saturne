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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_scheme_geometry.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for a face the weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Compute also the volume pefc attached to each edge of the face
 *
 * \param[in]       f          id of the face in the cell-wise numbering
 * \param[in]       pfq        primal face quantities
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       hf_coef    coefficient related to the height of p_{f,c}
 * \param[in]       f_coef     coefficient related to the area of f
 * \param[in, out]  wvf        pointer to an array storing the weight/vertex
 * \param[in, out]  pefc_vol   pointer to an array storing the volume of pefc

 */
/*----------------------------------------------------------------------------*/

inline static void
_get_wvf_pefcvol(short int                 f,
                 const cs_quant_t          pfq,
                 const cs_cell_mesh_t     *cm,
                 const double              hf_coef,
                 const double              f_coef,
                 cs_real_t                *wvf,
                 cs_real_t                *pefc_vol)
{
  double  xef_len;
  cs_real_3_t  xef_un, cp;

  assert(hf_coef > 0);

  /* Reset weights */
  for (int ii = 0; ii < cm->n_vc; ii++) wvf[ii] = 0;

  /* Compute a weight for each vertex of the current face */
  for (int ii = 0, i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++, ii++) {

    const short int  e = cm->f2e_ids[i];
    const cs_quant_t  peq = cm->edge[e];

    cs_math_3_length_unitv(peq.center, pfq.center, &xef_len, xef_un);
    cs_math_3_cross_product(xef_un, peq.unitv, cp);

    const double  tef = 0.5 * xef_len * peq.meas * cs_math_3_norm(cp);
    const double  ef_contrib = tef * f_coef; // tef = s(v1,e,f) + s(v2, e, f)

    pefc_vol[ii] = tef * hf_coef;
    wvf[cm->e2v_ids[2*e]] += ef_contrib;     // for v1
    wvf[cm->e2v_ids[2*e+1]] += ef_contrib;   // for v2

  } /* End of loop on face edges */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
  assert(fabs(hv) > cs_math_get_machine_epsilon()); /* Sanity check */

  const double  ohv1 = 1/hv;
  for (int k = 0; k < 3; k++)
    grd_v1[k] = unormal[k] * ohv1;

  /* Gradient for v2
     Normal direction to the plane in opposition to v2
     Height from this plane to the vertex v2 */
  cs_math_3_cross_product(uvc[v1], deq.unitv, unormal);
  hv = lvc[v2] * _dp3(uvc[v2], unormal);
  assert(fabs(hv) > cs_math_get_machine_epsilon()); /* Sanity check */

  const double  ohv2 = 1/hv;
  for (int k = 0; k < 3; k++)
    grd_v2[k] = unormal[k] * ohv2;
}

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
                    cs_real_t                *grd_c)
{
  double  xef_len;
  cs_real_3_t  xef_un, un;

  const cs_quant_t  pfq = cm->face[f];
  const cs_nvec3_t  deq = cm->dedge[f];

  /* Compute the area of each triangle of base e and apex f */
  for (int i = cm->f2e_idx[f], ii = 0; i < cm->f2e_idx[f+1]; i++, ii++) {

    const short int  e = cm->f2e_ids[i];
    const cs_quant_t  peq = cm->edge[e];

    cs_math_3_length_unitv(peq.center, pfq.center, &xef_len, xef_un);
    cs_math_3_cross_product(xef_un, peq.unitv, un);

    /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
    tef[ii] = 0.5 * xef_len*peq.meas * cs_math_3_norm(un);

  } /* End of loop on face edges */

  /* Compute the gradient of the Lagrange function related to xc
     which is constant inside p_{f,c} */
  const double  hf = _dp3(pfq.unitv, deq.unitv) * deq.meas;
  const cs_real_t  ohf = -cm->f_sgn[f]/hf;
  for (int k = 0; k < 3; k++)
    grd_c[k] = ohf * pfq.unitv[k];
}

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
                   cs_real_t                *pefc_vol)
{
  const cs_quant_t  pfq = cm->face[f];
  const cs_nvec3_t  deq = cm->dedge[f];
  const double  h_coef = cs_math_onethird * _dp3(pfq.unitv, deq.unitv)*deq.meas;
  const double  f_coef = 0.5/pfq.meas;

  assert(h_coef > 0);

  /* Compute geometric quantities */
  _get_wvf_pefcvol(f, pfq, cm, h_coef, f_coef, wvf, pefc_vol);
}

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
                   cs_real_t                *pefc_vol)
{
  const cs_quant_t  pfq = cm->face[f];
  const cs_nvec3_t  deq = cm->dedge[f];
  const double  h_coef = cs_math_onethird * _dp3(pfq.unitv, deq.unitv)*deq.meas;
  const double  f_coef = 0.5/pfq.meas;

  /* Compute geometric quantities */
  _get_wvf_pefcvol(f, pfq, cm, h_coef, f_coef, wvf, pefc_vol);

  return  h_coef * pfq.meas; // volume of p_{f,c}
}

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
cs_compute_fwbs_q2(short int                f,
                   const cs_cell_mesh_t    *cm,
                   cs_real_3_t              grd_c,
                   cs_real_t               *wvf,
                   cs_real_t               *pefc_vol)
{
  const cs_quant_t  pfq = cm->face[f];
  const cs_nvec3_t  deq = cm->dedge[f];
  const double  hf = _dp3(pfq.unitv, deq.unitv) * deq.meas;
  const double  f_coef = 0.5/pfq.meas;

  assert(hf > 0);

  /* Compute geometric quantities */
  _get_wvf_pefcvol(f, pfq, cm, cs_math_onethird * hf, f_coef, wvf, pefc_vol);

  /* Compute the gradient of the Lagrange function related to xc
     which is constant inside p_{f,c} */
  const cs_real_t  ohf = -cm->f_sgn[f]/hf;
  for (int k = 0; k < 3; k++)
    grd_c[k] = ohf * pfq.unitv[k];
}

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
                   cs_real_t                *pefc_vol)
{
  const cs_quant_t  pfq = cm->face[f];
  const cs_nvec3_t  deq = cm->dedge[f];
  const double  hf = _dp3(pfq.unitv, deq.unitv)*deq.meas;
  const double  h_coef = cs_math_onethird * hf;
  const double  f_coef = 0.5/pfq.meas;

  /* Compute geometric quantities */
  _get_wvf_pefcvol(f, pfq, cm, h_coef, f_coef, wvf, pefc_vol);

  /* Compute the gradient of the Lagrange function related to xc
     which is constant inside p_{f,c} */
  const cs_real_t  ohf = -cm->f_sgn[f]/hf;
  for (int k = 0; k < 3; k++)
    grd_c[k] = ohf * pfq.unitv[k];

  return  h_coef * pfq.meas; // volume of p_{f,c}
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
