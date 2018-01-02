#ifndef __CS_RECO_H__
#define __CS_RECO_H__

/*============================================================================
 * Routines to handle the reconstruction of fields
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"

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
 * \brief  Reconstruct at cell centers and face centers a vertex-based field
 *         Linear interpolation. If p_crec and/or p_frec are not allocated, this
 *         done in this subroutine.
 *
 *  \param[in]      connect  pointer to the connectivity struct.
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      dof      pointer to the field of vtx-based DoFs
 *  \param[in, out] p_crec   reconstructed values at cell centers
 *  \param[in, out] p_frec   reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_conf_vtx_dofs(const cs_cdo_connect_t      *connect,
                      const cs_cdo_quantities_t   *quant,
                      const double                *dof,
                      double                      *p_crec[],
                      double                      *p_frec[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at all cell centers from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   values of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pv_at_cell_centers(const cs_connect_index_t    *c2v,
                           const cs_cdo_quantities_t   *quant,
                           const double                *array,
                           cs_real_t                   *val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the cell center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      c_id     cell id
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pv_at_cell_center(cs_lnum_t                    c_id,
                          const cs_connect_index_t    *c2v,
                          const cs_cdo_quantities_t   *quant,
                          const double                *array,
                          cs_real_t                   *val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the face center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      fm     pointer to cs_face_mesh_t structure
 *  \param[in]      p_v    pointer to an array of values (local to this face)
 *  \param[in, out] p_f    value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_reco_pv_at_face_center(const cs_face_mesh_t    *fm,
                          const double            *p_v,
                          double                  *p_f)
{
  *p_f = 0.;

  if (p_v == NULL)
    return;

  const cs_quant_t  pfq = fm->face;

  for (short int e = 0; e < fm->n_ef; e++)
    *p_f += (p_v[fm->e2v_ids[2*e]] + p_v[fm->e2v_ids[2*e+1]]) * fm->tef[e];
  *p_f *= 0.5 / pfq.meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the face center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      f_id     face id (interior and border faces)
 *  \param[in]      connect  pointer to a cs_cdo_connect_t structure
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      pdi      pointer to the array of values
 *  \param[in, out] pdi_f    value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pf_from_pv(cs_lnum_t                     f_id,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant,
                   const double                 *pdi,
                   cs_real_t                    *pdi_f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector at the cell center from an array of
 *         values defined on dual faces lying inside each cell.
 *         This array is scanned thanks to the c2e connectivity.
 *
 *  \param[in]      c_id     cell id
 *  \param[in]      c2e      cell -> edges connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_at_cell_center(cs_lnum_t                    c_id,
                             const cs_connect_index_t    *c2e,
                             const cs_cdo_quantities_t   *quant,
                             const double                *array,
                             cs_real_3_t                  val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_c     value of the reconstructed vector in the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_cell(const cs_cell_mesh_t        *cm,
                      const cs_real_t             *array,
                      cs_real_3_t                  val_c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside pec which is a volume
 *         surrounding the edge e inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      e         local edge id
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_pec   value of the reconstruction in pec
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_pec(const cs_cell_mesh_t        *cm,
                     short int                    e,
                     const cs_real_t             *array,
                     cs_real_3_t                  val_pec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi          array of DoF values at vertices
 * \param[in]      length_xcxp  lenght of the segment [x_c, x_p]
 * \param[in]      unitv_xcxp   unitary vector pointed from x_c to x_p
 * \param[in, out] wbuf         pointer to a temporary buffer
 *
 * \return the value of the reconstructed potential at the evaluation point
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_reco_pv_inside_cell(const cs_cell_mesh_t    *cm,
                       const cs_real_t          pdi[],
                       const cs_real_t          length_xcxp,
                       const cs_real_t          unitv_xcxp[],
                       cs_real_t                wbuf[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct by a constant vector a field of edge-based DoFs
 *         in a volume surrounding an edge
 *
 *  \param[in]      cid     cell id
 *  \param[in]      e1_id   sub-volume related to this edge id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field in this sub-volume
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cost_edge_dof(cs_lnum_t                    cid,
                      cs_lnum_t                    e1_id,
                      const cs_connect_index_t    *c2e,
                      const cs_cdo_quantities_t   *quant,
                      const double                *dof,
                      double                       reco[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the cell center of the gradient of a field
 *         defined on primal vertices.
 *
 * \param[in]      c_id    cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant   pointer to the additional quantities struct.
 * \param[in]      pdi     pointer to the array of values
 * \param[in, out] val_xc  value of the reconstructed graident at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grd_cell_from_pv(cs_lnum_t                    c_id,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const double                *pdi,
                         cs_real_t                    val_xc[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at the cell center a field of edge-based DoFs
 *
 *  \param[in]      cid     cell id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dof(cs_lnum_t                   cid,
                      const cs_connect_index_t   *c2e,
                      const cs_cdo_quantities_t  *quant,
                      const double               *dof,
                      double                      reco[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at each cell center a field of edge-based DoFs
 *
 *  \param[in]      connect   pointer to the connectivity struct.
 *  \param[in]      quant     pointer to the additional quantities struct.
 *  \param[in]      dof       pointer to the field of edge-based DoFs
 *  \param[in, out] p_ccrec   pointer to the reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dofs(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *quant,
                       const double               *dof,
                       double                     *p_ccrec[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RECO_H__ */
