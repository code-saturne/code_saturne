#ifndef __CS_CDO_DIFFUSION_H__
#define __CS_CDO_DIFFUSION_H__

/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based  and vertex+cell schemes
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_hodge.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_cdo_diff_t  cs_cdo_diff_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure used to build the stiffness matrix
 *
 * \param[in] connect      pointer to a cs_cdo_connect_t structure
 * \param[in] space_scheme  scheme used for discretizing in space
 * \param[in] is_uniform   diffusion tensor is uniform ? (true or false)
 * \param[in] h_info       cs_param_hodge_t structure
 * \param[in] bc_enforce   type of boundary enforcement for Dirichlet values
 *
 * \return a pointer to a new allocated cs_cdo_diff_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_diff_t *
cs_cdo_diffusion_builder_init(const cs_cdo_connect_t       *connect,
                              cs_space_scheme_t             space_scheme,
                              bool                          is_uniform,
                              const cs_param_hodge_t        h_info,
                              const cs_param_bc_enforce_t   bc_enforce);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdo_diff_t structure
 *
 * \param[in, out ] diff   pointer to a cs_cdo_diff_t struc.
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

cs_cdo_diff_t *
cs_cdo_diffusion_builder_free(cs_cdo_diff_t   *diff);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the related Hodge builder structure
 *
 * \param[in]  diff   pointer to a cs_cdo_diff_t structure
 *
 * \return  a pointer to a cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_cdo_diffusion_get_hodge_builder(cs_cdo_diff_t   *diff);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get temporary buffers attached to a cs_cdo_diff_t structure
 *
 * \param[in]       diff     pointer to a cs_cdo_diff_t structure
 * \param[in, out]  tmp_vec  pointer to a buffer of cs_real_3_t
 * \param[in, out]  tmp_sca  pointer to a buffer of cs_real_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_get_tmp_buffers(const cs_cdo_diff_t   *diff,
                                 cs_real_3_t          **tmp_vec,
                                 cs_real_t            **tmp_sca);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a cell --> Dirichlet boundary faces connectivity
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      dir_face     pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] c2bcbf_idx   pointer to the index to build
 * \param[in, out] c2bcbf_ids   pointer to the list of ids to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_build_c2bcbf(const cs_cdo_connect_t    *connect,
                              const cs_cdo_bc_list_t    *dir_face,
                              cs_lnum_t                 *p_c2bcbf_idx[],
                              cs_lnum_t                 *p_c2bcbf_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) stiffness matrix
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      lm          cell-wise connectivity and quantitites
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdo_diffusion_build_local(const cs_cdo_quantities_t   *quant,
                             const cs_cell_mesh_t        *lm,
                             const cs_real_3_t           *tensor,
                             cs_cdo_diff_t               *diff);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) "normal trace gradient" matrix taking
 *          into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique (symmetrized or not)
 *
 * \param[in]       f_id      face id (a border face attached to a Dir. BC)
 * \param[in]       cm        pointer to a cs_cell_mesh_t struct.
 * \param[in]       matpty    3x3 matrix related to the diffusion property
 * \param[in, out]  diff      auxiliary structure used to build the diff. term
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_weak_bc(cs_lnum_t                    f_id,
                         const cs_cell_mesh_t        *cm,
                         const cs_real_t              matpty[3][3],
                         cs_cdo_diff_t               *diff,
                         cs_cdo_locsys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]       cm        pointer to a cs_face_mesh_t structure
 * \param[in]       dfaces    pointer to the dual faces related to cell edges
 * \param[in]       pty_tens  3x3 matrix related to the diffusion property
 * \param[in]       p_v       array of values attached to face vertices
 * \param[in]       p_c       value attached to the cell
 * \param[in, out]  diff      auxiliary structure dedicated to diffusion
 * \param[in, out]  c_flux    flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_cellwise_flux(const cs_cell_mesh_t      *cm,
                               const cs_dface_t          *dfaces,
                               const cs_real_t            pty_tens[3][3],
                               const double              *p_v,
                               const double               p_c,
                               cs_cdo_diff_t             *diff,
                               double                    *c_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across a face (based on a subdivision
 *          into tetrahedra of the volume p_{f,c})
 *
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       pty_tens  3x3 matrix related to the diffusion property
 * \param[in]       p_v       array of values attached to face vertices
 * \param[in]       p_f       value attached to the face
 * \param[in]       p_c       value attached to the cell
 * \param[in, out]  diff      auxiliary structure dedicated to diffusion
 *
 * \return the value of the diffusive flux across the current face
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_diffusion_face_flux(const cs_face_mesh_t      *fm,
                           const cs_real_t            pty_tens[3][3],
                           const double              *p_v,
                           const double               p_f,
                           const double               p_c,
                           cs_cdo_diff_t             *diff);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_DIFFUSION_H__ */
