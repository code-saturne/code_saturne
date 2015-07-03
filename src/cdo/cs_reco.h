#ifndef __CS_RECO_H__
#define __CS_RECO_H__

/*============================================================================
 * Routines to handle the reconstruction of fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include "cs_cdo_quantities.h"
#include "cs_cdo_connect.h"

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
 *  \param[in]    connect  pointer to the connectivity struct.
 *  \param[in]    quant    pointer to the additional quantities struct.
 *  \param[in]    dof      pointer to the field of vtx-based DoFs
 *  \param[inout] p_crec   reconstructed values at cell centers
 *  \param[inout] p_frec   reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_conf_vtx_dofs(const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      const double                 *dof,
                      double                       *p_crec[],
                      double                       *p_frec[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct by a constant vector a field of edge-based DoFs
 *         in a volume surrounding an edge
 *
 *  \param[in]    cid     cell id
 *  \param[in]    e1_id   sub-volume related to this edge id
 *  \param[in]    c2e     cell -> edges connectivity
 *  \param[in]    quant   pointer to the additional quantities struct.
 *  \param[in]    dof     pointer to the field of edge-based DoFs
 *  \param[inout] reco    value of the reconstrcuted field in this sub-volume
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dga_edge_dof(cs_lnum_t                    cid,
                     cs_lnum_t                    e1_id,
                     const cs_connect_index_t    *c2e,
                     const cs_cdo_quantities_t   *quant,
                     const double                *dof,
                     double                       reco[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at the cell center a field of edge-based DoFs
 *
 *  \param[in]    cid     cell id
 *  \param[in]    c2e     cell -> edges connectivity
 *  \param[in]    quant   pointer to the additional quantities struct.
 *  \param[in]    dof     pointer to the field of edge-based DoFs
 *  \param[inout] reco    value of the reconstrcuted field at cell center
 *
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
 *  \param[in]    connect   pointer to the connectivity struct.
 *  \param[in]    quant     pointer to the additional quantities struct.
 *  \param[in]    dof       pointer to the field of edge-based DoFs
 *  \param[inout] p_ccrec   pointer to the reconstructed values
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
