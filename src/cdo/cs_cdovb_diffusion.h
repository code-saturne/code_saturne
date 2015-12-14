#ifndef __CS_CDOVB_DIFFUSION_H__
#define __CS_CDOVB_DIFFUSION_H__

/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"
#include "cs_param.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_cdovb_diff_t  cs_cdovb_diff_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure used to build the stiffness matrix
 *
 * \param[in] connect      pointer to a cs_cdo_connect_t struct.
 * \param[in] is_uniform   diffusion tensor is uniform ? (true or false)
 * \param[in] h_info       cs_param_hodge_t struct.
 * \param[in] bc_enforce   type of boundary enforcement for Dirichlet values
 *
 * \return a pointer to a new allocated cs_cdovb_diffusion_builder_t struc.
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_diff_t *
cs_cdovb_diffusion_builder_init(const cs_cdo_connect_t       *connect,
                                bool                          is_uniform,
                                const cs_param_hodge_t        h_info,
                                const cs_param_bc_enforce_t   bc_enforce);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdovb_diff_t structure
 *
 * \param[in, out ] diff   pointer to a cs_cdovb_diff_t struc.
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_diff_t *
cs_cdovb_diffusion_builder_free(cs_cdovb_diff_t   *diff);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) stiffness matrix
 *
 * \param[in]      c_id        current cell id
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      vtag        pointer to a cs_cdovb_scaleq_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_diffusion_build_local(cs_lnum_t                    c_id,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_lnum_t             *vtag,
                               const cs_real_3_t           *tensor,
                               cs_cdovb_diff_t             *diff);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) "normal trace gradient" matrix
 *          This local matrix is used in Nitsche method to weakly penalized
 *          Dirichlet boundary conditions.
 *
 * \param[in]      c_id        cell id
 * \param[in]      f_id        face id (a border face attached to a Dir. BC)
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      matpty      3x3 matrix related to the diffusion property
 * \param[in]      eig_ratio   eigenvalue_max/eigenvalue_min
 * \param[in]      eig_max     eigenvalue with maximal value
 * \param[in, out] loc_v_ids   store local vertex ids
 * \param[in, out] v_coef      store local contribution on each border vertex
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local "normal trace gradient" matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_diffusion_ntrgrd_build(cs_lnum_t                    c_id,
                                cs_lnum_t                    f_id,
                                const cs_cdo_connect_t      *connect,
                                const cs_cdo_quantities_t   *quant,
                                const cs_real_t              matpty[3][3],
                                cs_real_t                    eig_ratio,
                                cs_real_t                    eig_max,
                                cs_lnum_t                   *loc_v_ids,
                                cs_real_t                   *v_coef,
                                cs_cdovb_diff_t             *diff);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_DIFFUSION_H__ */
