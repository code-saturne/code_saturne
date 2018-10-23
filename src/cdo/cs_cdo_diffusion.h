#ifndef __CS_CDO_DIFFUSION_H__
#define __CS_CDO_DIFFUSION_H__

/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based  and vertex+cell schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusion flux operator which corresponds to the normal
 *          trace operator for a given border face.
 *          Different algorithm can be used to reconstruct this flux.
 *
 * \param[in]      fm      pointer to a cs_face_mesh_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu     property tensor times the face normal
 * \param[in]      beta    value of the stabilization coef. or zero if not used
 * \param[in, out] cb      pointer to a cell builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op. i.e.
 *                         the flux operator
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_diffusion_flux_trace_t)(const cs_face_mesh_t     *fm,
                                const cs_cell_mesh_t     *cm,
                                const cs_real_3_t         mnu,
                                double                    beta,
                                cs_cell_builder_t        *cb,
                                cs_sdm_t                 *ntrgrd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Cellwise computation of the diffusive flux
 *
 * \param[in]      cm       pointer to a cs_face_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_cellwise_diffusion_flux_t)(const cs_cell_mesh_t     *cm,
                                   const cs_real_t          *pot,
                                   cs_cell_builder_t        *cb,
                                   cs_real_t                *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Enforce the Dirichlet BCs
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_diffusion_enforce_dir_t)(const cs_equation_param_t      *eqp,
                                 const cs_cell_mesh_t           *cm,
                                 cs_cdo_diffusion_flux_trace_t  *flux_op,
                                 cs_face_mesh_t                 *fm,
                                 cs_cell_builder_t              *cb,
                                 cs_cell_sys_t                  *csys);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique - Face-based version
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_diffusion_weak_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment - Face-based version
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_diffusion_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *          Specific to CDO-V+C schemes
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in]      beta      not useful here (prototype of function pointer)
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_diffusion_flux_op(const cs_face_mesh_t     *fm,
                            const cs_cell_mesh_t     *cm,
                            const cs_real_3_t         pty_nuf,
                            double                    beta,
                            cs_cell_builder_t        *cb,
                            cs_sdm_t                 *ntrgrd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in]      beta      not useful here (prototype of function pointer)
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_wbs_flux_op(const cs_face_mesh_t     *fm,
                               const cs_cell_mesh_t     *cm,
                               const cs_real_3_t         pty_nuf,
                               double                    beta,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *ntrgrd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          COST algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm      pointer to a cs_face_mesh_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu     property tensor times the face normal
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] cb      pointer to a cell builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op. i.e.
 *                         the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_cost_flux_op(const cs_face_mesh_t     *fm,
                                const cs_cell_mesh_t     *cm,
                                const cs_real_3_t         mnu,
                                double                    beta,
                                cs_cell_builder_t        *cb,
                                cs_sdm_t                 *ntrgrd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by keeping the DoFs related to
 *          Dirichlet BCs in the algebraic system (i.e. a weak enforcement)
 *          The corresponding DoFs are algebraically "removed" of the system
 *
 *          |      |     |     |      |     |     |  |     |          |
 *          | Aii  | Aid |     | Aii  |  0  |     |bi|     |bi-Aid.bd |
 *          |------------| --> |------------| and |--| --> |----------|
 *          |      |     |     |      |     |     |  |     |          |
 *          | Adi  | Add |     |  0   |  Id |     |bd|     |    xd    |
 *
 * where xd is the value of the Dirichlet BC
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_cdo_diffusion_flux_trace_t   *flux_op,
                                cs_face_mesh_t                  *fm,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by keeping the DoFs related to
 *          Dirichlet BCs in the algebraic system (i.e. a weak enforcement)
 *          The corresponding DoFs are algebraically "removed" of the system
 *          Block version.
 *
 *          |      |     |     |      |     |     |  |     |          |
 *          | Aii  | Aid |     | Aii  |  0  |     |bi|     |bi-Aid.xd |
 *          |------------| --> |------------| and |--| --> |----------|
 *          |      |     |     |      |     |     |  |     |          |
 *          | Adi  | Add |     |  0   |  Id |     |bd|     |    xd    |
 *
 * where xd is the value of the Dirichlet BC
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_cdo_diffusion_flux_trace_t   *flux_op,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_cdo_diffusion_flux_trace_t   *flux_op,
                                cs_face_mesh_t                  *fm,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value.
 *          Case of a cellwise system defined by block.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_cdo_diffusion_flux_trace_t   *flux_op,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_weak_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_wsym_dirichlet(const cs_equation_param_t       *eqp,
                                  const cs_cell_mesh_t            *cm,
                                  cs_cdo_diffusion_flux_trace_t   *flux_op,
                                  cs_face_mesh_t                  *fm,
                                  cs_cell_builder_t               *cb,
                                  cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          Use the COST algo. for computing the discrete Hodge op.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Hdg * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcost_get_dfbyc_flux(const cs_cell_mesh_t      *cm,
                                      const double              *pot,
                                      cs_cell_builder_t         *cb,
                                      double                    *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux inside a (primal) cell
 *          Use the COST algo. for computing the discrete Hodge op.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Hdg * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcost_get_pc_flux(const cs_cell_mesh_t      *cm,
                                   const double              *pot,
                                   cs_cell_builder_t         *cb,
                                   double                    *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_dfbyc_flux(const cs_cell_mesh_t   *cm,
                                    const cs_real_t        *pot,
                                    cs_cell_builder_t      *cb,
                                    cs_real_t              *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux inside a given primal cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux vector inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_pc_flux(const cs_cell_mesh_t   *cm,
                                 const cs_real_t        *pot,
                                 cs_cell_builder_t      *cb,
                                 cs_real_t              *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing a piecewise constant gradient from the degrees of
 *          freedom.
 *
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_tensor  property tensor times the face normal
 * \param[in]      pot_values   array of values of the potential (all the mesh)
 * \param[in]      f            face id in the cell mesh
 * \param[in]      t_eval       time at which one evaluates the advection field
 * \param[in, out] fluxes       values of the fluxes related to each vertex
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_face_p0_flux(const cs_cell_mesh_t     *cm,
                                const cs_real_3_t        *diff_tensor,
                                const cs_real_t          *pot_values,
                                short int                 f,
                                cs_real_t                 t_eval,
                                cs_real_t                *fluxes);

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
 * \param[in, out]  cb        auxiliary structure dedicated to diffusion
 *
 * \return the value of the diffusive flux across the current face
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_diffusion_face_wbs_flux(const cs_face_mesh_t      *fm,
                               const cs_real_t            pty_tens[3][3],
                               const double              *p_v,
                               const double               p_f,
                               const double               p_c,
                               cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing the degrees of freedom.
 *          This routine shares similarities with
 *          \ref cs_cdovb_diffusion_cost_flux_op
 *
 * \param[in]   f             face id in the cell mesh
 * \param[in]   cm            pointer to a cs_cell_mesh_t structure
 * \param[in]   diff_tensor   property tensor times the face normal
 * \param[in]   pot_values    array of values of the potential (all the mesh)
 * \param[in]   beta          value of the stabilization coef. related to reco.
 * \param[in, out]  cb        auxiliary structure dedicated to diffusion
 *
 * \return the diffusive flux across the face f
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdovb_diffusion_face_cost_flux(short int                 f,
                                  const cs_cell_mesh_t     *cm,
                                  const cs_real_3_t        *diff_tensor,
                                  const cs_real_t          *pot_values,
                                  double                    beta,
                                  cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_DIFFUSION_H__ */
