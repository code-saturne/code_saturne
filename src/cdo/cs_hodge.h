#ifndef __CS_HODGE_H__
#define __CS_HODGE_H__

/*============================================================================
 * Manage discrete Hodge operators and closely related operators
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_property.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local operator for a given cell (stored in cb->hdg for a
 *          discrete Hodge operator or in cb->loc for a stiffness matrix)
 *
 * \param[in]      hodgep    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_hodge_t)(const cs_param_hodge_t    hodgep,
             const cs_cell_mesh_t     *cm,
             cs_cell_builder_t        *cb);

/*============================================================================
 * Static inline public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check if two sets of parameters related to how build a discrete
 *          Hodge operator are similar.
 *
 * \param[in]  h1_info     pointer to a first cs_param_hodge_t structure
 * \param[in]  h2_info     pointer to a second cs_param_hodge_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_hodge_param_is_similar(const cs_param_hodge_t    h1_info,
                          const cs_param_hodge_t    h2_info)
{
  if (h1_info.type != h2_info.type)
    return false;
  if (h1_info.algo != h2_info.algo)
    return false;
  if (h1_info.algo == CS_PARAM_HODGE_ALGO_COST ||
      h1_info.algo == CS_PARAM_HODGE_ALGO_BUBBLE) {
    if (fabs(h1_info.coef - h2_info.coef) > 0)
      return false;
    else
      return true;
  }
  else
    return true;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc
 *          Case of CDO face-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_cost_get_stiffness(const cs_param_hodge_t    hodgep,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          with the usage of bubble stabilization.
 *          The computed matrix is stored in cb->loc
 *          Case of CDO face-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_bubble_get_stiffness(const cs_param_hodge_t    hodgep,
                                 const cs_cell_mesh_t     *cm,
                                 cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc
 *          Case of CDO face-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_voro_get_stiffness(const cs_param_hodge_t    hodgep,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_cost_get_stiffness(const cs_param_hodge_t    hodgep,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          Case of CDO vertex-based schemes and isotropic property
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_cost_get_iso_stiffness(const cs_param_hodge_t    hodgep,
                                   const cs_cell_mesh_t     *cm,
                                   cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_cost_get_aniso_stiffness(const cs_param_hodge_t    hodgep,
                                     const cs_cell_mesh_t     *cm,
                                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          Case of CDO vertex-based schemes and isotropic material property
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_bubble_get_iso_stiffness(const cs_param_hodge_t    hodgep,
                                     const cs_cell_mesh_t     *cm,
                                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          Case of CDO vertex-based schemes and anisotropic material property
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_bubble_get_aniso_stiffness(const cs_param_hodge_t    hodgep,
                                       const cs_cell_mesh_t     *cm,
                                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          Case of anisotropic material property
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_ocs2_get_aniso_stiffness(const cs_param_hodge_t    hodgep,
                                     const cs_cell_mesh_t     *cm,
                                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_voro_get_stiffness(const cs_param_hodge_t    hodgep,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS) algo.
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_wbs_get_stiffness(const cs_param_hodge_t    hodgep,
                              const cs_cell_mesh_t     *cm,
                              cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS) algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      hodgep     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vcb_get_stiffness(const cs_param_hodge_t    hodgep,
                           const cs_cell_mesh_t     *cm,
                           cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator on a given cell which is equivalent of
 *         a mass matrix. It relies on a CO+ST algo. and is specific to CDO-Fb
 *         schemes.
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_get_mass(const cs_param_hodge_t    hodgep,
                     const cs_cell_mesh_t     *cm,
                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the WBS algo.
 *          This function is specific for vertex+cell-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vcb_wbs_get(const cs_param_hodge_t    hodgep,
                     const cs_cell_mesh_t     *cm,
                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using WBS algo.
 *          Hodge op. from primal vertices to dual cells.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vpcd_wbs_get(const cs_param_hodge_t    hodgep,
                      const cs_cell_mesh_t     *cm,
                      cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal vertices to dual cells.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vpcd_voro_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_voro_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_cost_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          with a bubble stabilization.
 *          Hodge op. from primal edges to dual faces. This function is
 *          specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_bubble_get(const cs_param_hodge_t    hodgep,
                         const cs_cell_mesh_t     *cm,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_ocs2_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal faces to dual edges.
 *          This function is related to cell-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fped_voro_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal faces to dual edges.
 *          This function is related to cell-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fped_cost_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from primal faces to dual edges.
 *          This function is related to cell-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fped_bubble_get(const cs_param_hodge_t    hodgep,
                         const cs_cell_mesh_t     *cm,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_voro_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_cost_get(const cs_param_hodge_t    hodgep,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_bubble_get(const cs_param_hodge_t    hodgep,
                         const cs_cell_mesh_t     *cm,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_cost_get_opt(const cs_param_hodge_t    hodgep,
                           const cs_cell_mesh_t     *cm,
                           cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator from dual edges to primal faces for a
 *          given cell using the COST algo. with a bubble stabilization.
 *          This function is related to face-based schemes
 *
 * \param[in]      hodgep      pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_bubble_get(const cs_param_hodge_t    hodgep,
                         const cs_cell_mesh_t     *cm,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator and multiple it with
 *          a vector
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      hodgep      cs_param_hodge_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      in_vals   vector to multiply with the discrete Hodge op.
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] result    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_matvec(const cs_cdo_connect_t       *connect,
                const cs_cdo_quantities_t    *quant,
                const cs_param_hodge_t        hodgep,
                const cs_property_t          *pty,
                const cs_real_t               in_vals[],
                cs_real_t                     t_eval,
                cs_real_t                     result[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator in order to define
 *          a circulation array from a flux array
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      hodgep      cs_param_hodge_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      flux      vector to multiply with the discrete Hodge op.
 * \param[in, out] circul    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_circulation_from_flux(const cs_cdo_connect_t       *connect,
                               const cs_cdo_quantities_t    *quant,
                               cs_real_t                     t_eval,
                               const cs_param_hodge_t        hodgep,
                               const cs_property_t          *pty,
                               const cs_real_t               flux[],
                               cs_real_t                     circul[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the hodge operator related to a face (i.e. a mass matrix
 *          with unity property) using a Whitney Barycentric Subdivision (WBS)
 *          algorithm
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] hf        pointer to a cs_sdm_t structure to define
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_compute_wbs_surfacic(const cs_face_mesh_t    *fm,
                              cs_sdm_t                *hf);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HODGE_H__ */
