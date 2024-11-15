#ifndef __CS_RECO_H__
#define __CS_RECO_H__

/*============================================================================
 * Functions to handle the reconstruction of fields
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply 1/|dual_vol| to a synchronized array of DoF vertices
 *         Parallel synchronization is done inside this function:
 *         1) A parallel sum reduction is done on the vtx_values
 *         2) Apply 1/|dual_vol| to each entry
 *
 *  \param[in]      connect    pointer to additional connectivities for CDO
 *  \param[in]      quant      pointer to additional quantities for CDO
 *  \param[in]      stride     number of entries for each vertex
 *  \param[in]      interlace  interlaced array (useful if the stride > 1)
 *  \param[in, out] array      array of DoF values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dual_vol_weight_reduction(const cs_cdo_connect_t       *connect,
                                  const cs_cdo_quantities_t    *quant,
                                  int                           stride,
                                  bool                          interlace,
                                  cs_real_t                    *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]   c2f           cell -> faces connectivity
 * \param[in]   cdoq          pointer to the additional quantities struct.
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
 * \param[in]  c_id          id of the cell to handle
 * \param[in]  c2f           cell -> faces connectivity
 * \param[in]  cdoq          pointer to the additional quantities struct.
 * \param[in]  face_dofs     array of DoF values at faces
 * \param[in]  local_input   true means that face_dofs is of size n_cell_faces
 * \param[out] cell_reco     vector-valued reconstruction inside cells. This
 *                           quantity should have been allocated before calling
 *                           this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vector_by_face_dofs(cs_lnum_t                   c_id,
                                 const cs_adjacency_t       *c2f,
                                 const cs_cdo_quantities_t  *cdoq,
                                 const cs_real_t             face_dofs[],
                                 bool                        local_input,
                                 cs_real_t                  *cell_reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the face
 *        DoFs (interior and boundary). Scalar-valued face DoFs are related to
 *        the normal flux across faces.
 *
 * \param[in]  c2f          cell -> faces connectivity
 * \param[in]  cdoq         pointer to the additional quantities struct.
 * \param[in]  face_dofs    array of DoF values at faces
 * \param[out] cell_reco    vector-valued reconstruction inside cells. This
 *                          quantity should have been allocated before calling
 *                          this function
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cell_vectors_by_face_dofs(const cs_adjacency_t       *c2f,
                                  const cs_cdo_quantities_t  *cdoq,
                                  const cs_real_t             face_dofs[],
                                  cs_real_t                  *cell_reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        at primal vertices.
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or null
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2c(cs_lnum_t                    n_cells,
                   const cs_lnum_t             *cell_ids,
                   const cs_adjacency_t        *c2v,
                   const cs_cdo_quantities_t   *cdoq,
                   const double                *array,
                   bool                         dense_ouput,
                   cs_real_t                   *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (v, c) --> array which can be scanned by the c2v
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or null
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (v,c)
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_vbyc2c(cs_lnum_t                    n_cells,
                      const cs_lnum_t             *cell_ids,
                      const cs_adjacency_t        *c2v,
                      const cs_cdo_quantities_t   *cdoq,
                      const double                *array,
                      bool                         dense_ouput,
                      cs_real_t                   *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (e, c) --> array which can be scanned by the c2e
 *        adjacency).
 *        Case of scalar-valued array.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or null
 * \param[in]      c2e           cell -> edges connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (e,c)
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_ebyc2c(cs_lnum_t                    n_cells,
                      const cs_lnum_t             *cell_ids,
                      const cs_adjacency_t        *c2e,
                      const cs_cdo_quantities_t   *cdoq,
                      const double                *array,
                      bool                         dense_ouput,
                      cs_real_t                   *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at all cell centers from an array of values
 *        defined on primal vertices.
 *        Case of vector-valued fields.
 *
 * \param[in]      n_cells       number of selected cells
 * \param[in]      cell_ids      list of cell ids or null
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in]      dense_ouput   apply cell_ids on the reco array
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vector_v2c(cs_lnum_t                    n_cells,
                   const cs_lnum_t             *cell_ids,
                   const cs_adjacency_t        *c2v,
                   const cs_cdo_quantities_t   *cdoq,
                   const double                *array,
                   bool                         dense_ouput,
                   cs_real_t                   *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at the face center from an array of values
 *        defined on primal vertices.
 *        Case of scalar-valued arrays.
 *
 * \param[in]      n_faces      number of faces
 * \param[in]      face_ids     list of face ids (interior/border faces) or null
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq         pointer to the additional quantities struct.
 * \param[in]      array        pointer to the array of values at vertices
 * \param[in]      dense_ouput  apply cell_ids on the reco array
 * \param[in, out] reco         value of the reconstruction at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2f(cs_lnum_t                     n_faces,
                   const cs_lnum_t              *face_ids,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq,
                   const double                 *array,
                   bool                          dense_ouput,
                   cs_real_t                    *reco);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct at face centers by a cell-based field
 *        weighted average.
 *
 *  \param[in]      connect   pointer to additional connectivities for CDO
 *  \param[in]      cdoq      pointer to additional quantities for CDO
 *  \param[in]      p_c       dofs at cell centers
 *  \param[in, out] p_reco_f  reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_c2f(const cs_cdo_connect_t     *connect,
                   const cs_cdo_quantities_t  *cdoq,
                   const double               *p_c,
                   double                     *p_reco_f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct at cell centers and face centers a vertex-based field
 *        Linear interpolation. If p_reco_c and/or p_reco_f are not allocated,
 *        this is done in this subroutine.
 *        Case of scalar-valued arrays.
 *
 *  \param[in]      connect   pointer to additional connectivities for CDO
 *  \param[in]      cdoq      pointer to additional quantities for CDO
 *  \param[in]      dof       pointer to the field of vtx-based DoFs
 *  \param[in, out] p_reco_c  reconstructed values at cell centers
 *  \param[in, out] p_reco_f  reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scalar_v2c_v2f(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *cdoq,
                       const double               *dof,
                       double                     *p_reco_c[],
                       double                     *p_reco_f[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a scalar-valued array at vertices from a scalar-valued
 *         array at cells.
 *
 * \param[in]      connect   pointer to additional connectivities for CDO
 * \param[in]      quant     pointer to the additional quantities for CDO
 * \param[in]      cell_val  array of scalar-valued values at cells
 * \param[in, out] vtx_val   array of scalar-valued values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_scal_pv_from_pc(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_real_t             *cell_val,
                        cs_real_t                   *vtx_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a vector-valued array at vertices from a vector-valued
 *         array at cells.
 *
 * \param[in]      connect   pointer to additional connectivities for CDO
 * \param[in]      quant     pointer to the additional quantities for CDO
 * \param[in]      cell_val  array of vector-valued values at cells
 * \param[in, out] vtx_val   array of vector-valued values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_from_pc(const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_real_t             *cell_val,
                        cs_real_t                   *vtx_val);

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
 *                          cells (allocated if set to null)
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
 * \brief Reconstruct the constant gradient vector in a cell (the mean value)
 *        from the value at mesh vertices.
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      pdi      pointer to the array of values
 * \param[in, out] grdc     value of the reconstructed gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grad_cell_from_pv(cs_lnum_t                    c_id,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *pdi,
                          cs_real_t                    grdc[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the constant gradient vector in a cell (the mean value)
 *        from the value at mesh vertices. Case of two scalar fields.
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      p1di     pointer to the array of values
 * \param[in]      p2di     pointer to the array of values
 * \param[in, out] grd1c    value of the reconstructed gradient for p1
 * \param[in, out] grd2c    value of the reconstructed gradient for p2
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_2grad_cell_from_pv(cs_lnum_t                    c_id,
                           const cs_cdo_connect_t      *connect,
                           const cs_cdo_quantities_t   *quant,
                           const cs_real_t             *p1di,
                           const cs_real_t             *p2di,
                           cs_real_t                    grd1c[],
                           cs_real_t                    grd2c[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the vector-valued quantity inside each cell from the
 *        given flux array. This array is stored in the same order as cm->f_ids
 *        Scalar-valued face DoFs are related to the normal flux across primal
 *        faces.
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes         array of normal fluxes on primal faces
 * \param[out]     cell_reco      vector-valued reconstruction inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_vect_from_flux(const cs_cell_mesh_t    *cm,
                               const cs_real_t         *fluxes,
                               cs_real_t               *cell_reco);

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

/*============================================================================
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        at primal vertices.
 *        Case of scalar-valued array with all cells selected.
 *
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_reco_scalar_v2c_full(const cs_adjacency_t        *c2v,
                        const cs_cdo_quantities_t   *cdoq,
                        const double                *array,
                        cs_real_t                   *reco)
{
  cs_reco_scalar_v2c(cdoq->n_cells, NULL, c2v, cdoq, array, false, reco);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        at primal vertices.
 *        Case of vector-valued array with all cells selected.
 *
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at vertices
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_reco_vector_v2c_full(const cs_adjacency_t        *c2v,
                        const cs_cdo_quantities_t   *cdoq,
                        const double                *array,
                        cs_real_t                   *reco)
{
  cs_reco_vector_v2c(cdoq->n_cells, NULL, c2v, cdoq, array, false, reco);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (v, c) --> array which can be scanned by the c2v
 *        adjacency).
 *        Case of scalar-valued array with all cells selected.
 *
 * \param[in]      c2v           cell -> vertices connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (v,c)
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_reco_scalar_vbyc2c_full(const cs_adjacency_t        *c2v,
                           const cs_cdo_quantities_t   *cdoq,
                           const double                *array,
                           cs_real_t                   *reco)
{
  cs_reco_scalar_vbyc2c(cdoq->n_cells, NULL, c2v, cdoq, array, false, reco);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct the value at cell center from an array of values defined
 *        for each couple (e, c) --> array which can be scanned by the c2e
 *        adjacency).
 *        Case of scalar-valued array with all cells selected.
 *
 * \param[in]      c2e           cell -> edges connectivity
 * \param[in]      cdoq          pointer to the additional quantities struct.
 * \param[in]      array         pointer to the array of values at (e,c)
 * \param[in, out] reco          reconstructed values at the cell centers
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_reco_scalar_ebyc2c_full(const cs_adjacency_t        *c2e,
                           const cs_cdo_quantities_t   *cdoq,
                           const double                *array,
                           cs_real_t                   *reco)
{
  cs_reco_scalar_ebyc2c(cdoq->n_cells, NULL, c2e, cdoq, array, false, reco);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RECO_H__ */
