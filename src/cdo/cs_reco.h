#ifndef __CS_RECO_H__
#define __CS_RECO_H__

/*============================================================================
 * Routines to handle the reconstruction of fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "assert.h"

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_flag.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the face center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      fm     pointer to cs_face_mesh_t structure
 *  \param[in]      p_v    pointer to an array of values (local to this face)
 *
 * \return the value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_fw_scalar_pv_at_face_center(const cs_face_mesh_t    *fm,
                                    const cs_real_t         *p_v)
{
  cs_real_t  p_f = 0.;

  if (p_v == NULL)
    return p_f;

  const cs_quant_t  pfq = fm->face;

  for (short int e = 0; e < fm->n_ef; e++)
    p_f += (p_v[fm->e2v_ids[2*e]] + p_v[fm->e2v_ids[2*e+1]]) * fm->tef[e];
  p_f *= 0.5 / pfq.meas;

  return p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the face center from an array of values
 *         defined on primal vertices.
 *
 * \param[in]  f      id of the face in the cellwise numbering
 * \param[in]  cm     pointer to cs_cell_mesh_t structure
 * \param[in]  p_v    pointer to an array of values (local to the cell)
 *
 * \return the value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_pv_at_face_center(const short int          f,
                                    const cs_cell_mesh_t    *cm,
                                    const cs_real_t         *p_v)
{
  cs_real_t  p_f = 0.;

  if (p_v == NULL)
    return p_f;

  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ |
                       CS_FLAG_COMP_FE));

  for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {
    const short int  *v = cm->e2v_ids + 2*cm->f2e_ids[ie];
    p_f += (p_v[v[0]] + p_v[v[1]]) * cm->tef[ie];
  }
  p_f *= 0.5 / cm->face[f].meas;

  return p_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at the cell center from
 *         an array of values defined on primal vertices.
 *         Algorithm based on the cs_cell_mesh_t structure.
 *
 *  \param[in]      cm       pointer to a cs_cell_mesh_t structure
 *  \param[in]      array    pointer to the array of values (size: n_vc)
 *
 *  \return the value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_reco_cw_scalar_pv_at_cell_center(const cs_cell_mesh_t     *cm,
                                    const cs_real_t          *p_v)
{
  cs_real_t  p_c = 0.;

  if (p_v == NULL || cm == NULL)
    return p_c;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ)); /* Sanity check */

  /* Reconstruct the value at the cell center */
  for (short int v = 0; v < cm->n_vc; v++)
    p_c += cm->wvc[v] * p_v[v];

  return p_c;
}

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
cs_reco_pv_at_cell_centers(const cs_adjacency_t        *c2v,
                           const cs_cdo_quantities_t   *quant,
                           const double                *array,
                           cs_real_t                   *val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at all cell centers from an array of values
 *         defined on primal vertices.
 *         Case of vector-valued fields.
 *
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   values of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_at_cell_centers(const cs_adjacency_t        *c2v,
                                const cs_cdo_quantities_t   *quant,
                                const double                *array,
                                cs_real_t                   *val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]   c2f           cell -> faces connectivity
 * \param[in]   quant         pointer to the additional quantities struct.
 * \param[in]   i_face_vals   array of DoF values for interior faces
 * \param[in]   b_face_vals   array of DoF values for border faces
 * \param[out]  cell_reco     vector-valued reconstruction inside cells
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vectors_by_ib_face_dofs(const cs_adjacency_t       *c2f,
                                     const cs_cdo_quantities_t  *cdoq,
                                     const cs_real_t             i_face_vals[],
                                     const cs_real_t             b_face_vals[],
                                     cs_real_t                  *cell_reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside a cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]   c_id          id of the cell to handle
 * \param[in]   c2f           cell -> faces connectivity
 * \param[in]   quant         pointer to the additional quantities struct.
 * \param[in]   face_dofs     array of DoF values at faces
 * \param[out]  cell_reco     vector-valued reconstruction inside cells. This
 *                            quantity should have been allocated before calling
 *                            this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vector_by_face_dofs(cs_lnum_t                   c_id,
                                 const cs_adjacency_t       *c2f,
                                 const cs_cdo_quantities_t  *cdoq,
                                 const cs_real_t             face_dofs[],
                                 cs_real_t                  *cell_reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]   c2f           cell -> faces connectivity
 * \param[in]   quant         pointer to the additional quantities struct.
 * \param[in]   face_dofs     array of DoF values at faces
 * \param[out]  cell_reco     vector-valued reconstruction inside cells. This
 *                            quantity should have been allocated before calling
 *                            this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vectors_by_face_dofs(const cs_adjacency_t       *c2f,
                                  const cs_cdo_quantities_t  *cdoq,
                                  const cs_real_t             face_dofs[],
                                  cs_real_t                  *cell_reco);

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
                          const cs_adjacency_t        *c2v,
                          const cs_cdo_quantities_t   *quant,
                          const double                *array,
                          cs_real_t                   *val_xc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a vector-valued array at vertices from a vector-valued
 *         array at cells.
 *
 *  \param[in]      c2v       cell -> vertices connectivity
 *  \param[in]      quant     pointer to the additional quantities struct.
 *  \param[in]      val       pointer to the array of values
 *  \param[in, out] reco_val  values of the reconstruction at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_from_pc(const cs_adjacency_t        *c2v,
                        const cs_cdo_quantities_t   *quant,
                        const double                *val,
                        cs_real_t                   *reco_val);

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
                             const cs_adjacency_t        *c2e,
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
 * \brief  Reconstruct at the cell center a field of edge-based DoFs
 *
 *  \param[in]      c_id    cell id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dof(cs_lnum_t                   c_id,
                      const cs_adjacency_t       *c2e,
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
/*!
 * \brief  Reconstruct a cell-wise constant curl from the knowledge of the
 *         circulation at primal edges
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      circ     pointer to the array of circulations at edges
 * \param[in, out] p_curl   pointer to value of the reconstructed curl inside
 *                          cells (allocated if set to NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_curl_by_edge_dofs(const cs_cdo_connect_t        *connect,
                               const cs_cdo_quantities_t     *quant,
                               const cs_real_t               *circ,
                               cs_real_t                    **p_curl);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the mean-value of the gradient field with DoFs arising
 *         from a face-based scheme (values at face center and cell center)
 *         The reconstruction only deals with the consistent part so that there
 *         is no distinction betwwen Fb schemes
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      p_c      pointer to the array of values in cells
 * \param[in]      p_f      pointer to the array of values on faces
 * \param[in, out] grd_c    value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_cell_from_fb_dofs(cs_lnum_t                    c_id,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_real_t             *p_c,
                               const cs_real_t             *p_f,
                               cs_real_t                    grd_c[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the mean-value of the tensor gradient field with DoFs
 *         arising from a face-based scheme (vector-valued at face center and
 *         cell center) The reconstruction only deals with the consistent part
 *         so that there is no distinction between Fb schemes
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      u_c      pointer to the array of values in cells
 * \param[in]      u_f      pointer to the array of values on faces
 * \param[in, out] grd_c    value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_33_cell_from_fb_dofs(cs_lnum_t                    c_id,
                                  const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  const cs_real_t             *u_c,
                                  const cs_real_t             *u_f,
                                  cs_real_t                    grd_c[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the cell center of the gradient of a field
 *         defined on primal vertices.
 *
 * \param[in]      c_id    cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant   pointer to the additional quantities struct.
 * \param[in]      pdi     pointer to the array of values
 * \param[in, out] val_xc  value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_cell_from_pv(cs_lnum_t                    c_id,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *pdi,
                          cs_real_t                    val_xc[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      i_face_vals    array of DoF values for interior faces
 * \param[in]      b_face_vals    array of DoF values for border faces
 * \param[out]     cell_reco      vector-valued reconstruction inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_vect_from_face_dofs(const cs_cell_mesh_t    *cm,
                                    const cs_real_t          i_face_vals[],
                                    const cs_real_t          b_face_vals[],
                                    cs_real_t               *cell_reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi            array of DoF values at vertices
 * \param[out]     cell_gradient  gradient inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_grad_from_scalar_pv(const cs_cell_mesh_t    *cm,
                                    const cs_real_t          pdi[],
                                    cs_real_t               *cell_gradient);

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
cs_reco_cw_scalar_pv_inside_cell(const cs_cell_mesh_t    *cm,
                                 const cs_real_t          pdi[],
                                 const cs_real_t          length_xcxp,
                                 const cs_real_t          unitv_xcxp[],
                                 cs_real_t                wbuf[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the weighted (by volume) gradient inside a given primal
 *          cell for the related vertices.
 *          Use the WBS algo. for approximating the gradient.
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] vgrd     gradient at vertices inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_vgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *vgrd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the mean value of a gradient inside a given primal cell.
 *          Use the WBS algo. for approximating the gradient.
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] cgrd     mean-value of the cell gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *cgrd);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RECO_H__ */
