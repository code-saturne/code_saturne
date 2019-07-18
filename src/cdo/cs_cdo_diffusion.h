#ifndef __CS_CDO_DIFFUSION_H__
#define __CS_CDO_DIFFUSION_H__

/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based  and vertex+cell schemes
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
 * \brief   Cellwise computation of the diffusive flux
 *
 * \param[in]      cm       pointer to a cs_face_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_diffusion_cw_flux_t)(const cs_cell_mesh_t     *cm,
                             const cs_real_t          *pot,
                             cs_cell_builder_t        *cb,
                             cs_real_t                *flx);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
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
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys);

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
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
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
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique.
 *          Case of scalar-valued CDO Face-based schemes
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique.
 *          Case of vector-valued CDO Face-based schemes
 *          The idea is to compute the scalar version and dispatch it three
 *          times, one for each Cartesian components
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment - Face-based version
 *          Case of scalar-valued CDO Face-based schemes
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *          The idea is to compute the scalar version and dispatch it three
 *          times, one for each Cartesian components
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account sliding BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_wsym_sliding(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Robin BCs.
 *          Case of scalar-valued CDO-Vb schemes with a CO+ST algorithm.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_cost_robin(const cs_equation_param_t      *eqp,
                                const cs_cell_mesh_t           *cm,
                                cs_face_mesh_t                 *fm,
                                cs_cell_builder_t              *cb,
                                cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Robin BCs.
 *          Case of scalar-valued CDO-Vb schemes with a WBS algorithm.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_robin(const cs_equation_param_t      *eqp,
                               const cs_cell_mesh_t           *cm,
                               cs_face_mesh_t                 *fm,
                               cs_cell_builder_t              *cb,
                               cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account generic BCs by a weak enforcement using Nitsche
 *          technique. According to the settings one can apply Neumann BCs if
 *          alpha = 0, Dirichlet BCs if alpha >> 1 or Robin BCs
 *          Case of scalar-valued CDO-Vb schemes with a CO+ST algorithm.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_cost_generic(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of scalar-valued CDO-Vb schemes with an orthogonal
 *          splitting between the consistency/stabilization parts (OCS)
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_ocs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. A Dirichlet is set for the three components of the
 *          vector. Case of vector-valued CDO-Vb schemes with an orthogonal
 *          splitting between the consistency/stabilization parts (OCS)
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vvb_ocs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a sliding BCs.
 *          Case of vector-valued CDO-Vb schemes with a OCS algorithm.
 *          Orthogonal splitting betwwen Consistency/Stabilization parts.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vvb_ocs_sliding(const cs_equation_param_t      *eqp,
                                 const cs_cell_mesh_t           *cm,
                                 cs_face_mesh_t                 *fm,
                                 cs_cell_builder_t              *cb,
                                 cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-Vb schemes with a
 *          COST/Bubble or Voronoi algorithm. One assumes an Orthogonal
 *          splitting between Consistency/Stabilization parts (OCS).
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_ocs_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of CDO-Vb schemes with a WBS algorithm.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-Vb schemes with a
 *          WBS algorithm
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of CDO-VCb schemes with a WBS algorithm.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-VCb schemes with
 *          a WBS algorithm
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell.
 *          Use the same consistent approximation as in the discrete Hodge op.
 *          for this computation.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Consistent(Hdg) * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_get_dfbyc_flux(const cs_cell_mesh_t      *cm,
                                    const double              *pot,
                                    cs_cell_builder_t         *cb,
                                    double                    *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the constant approximation of the diffusive flux inside a
 *          (primal) cell. Use the same consistent approximation as in the
 *          discrete Hodge op. for this computation. This function is dedicated
 *          to vertex-based schemes. Flux = -Hdg * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_get_cell_flux(const cs_cell_mesh_t      *cm,
                                   const double              *pot,
                                   cs_cell_builder_t         *cb,
                                   double                    *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the normal flux for a face assuming only the knowledge of
 *         the potential at cell vertices. Valid for algorithm relying on a
 *         spliting of the consistency/stabilization part as in OCS (also
 *         called CO+ST) or Bubble algorithm. This is used for reconstructing
 *         the normal flux from the degrees of freedom. The contribution for
 *         each vertex of the face is then computed.
 *
 * \param[in]      f       face id in the cell mesh
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     array of values of the potential (all the mesh)
 * \param[in, out] cb      auxiliary structure dedicated to diffusion
 * \param[in, out] flux    array of values to set (size: n_vc)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_vbyf_flux(short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux);

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
cs_cdo_diffusion_wbs_get_cell_flux(const cs_cell_mesh_t   *cm,
                                   const cs_real_t        *pot,
                                   cs_cell_builder_t      *cb,
                                   cs_real_t              *flx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices (and at cell center).
 *          WBS algorithm is used for reconstructing the normal flux from the
 *          degrees of freedom.
 *
 * \param[in]  f              face id in the cell mesh
 * \param[in]  eqp            pointer to a cs_equation_param_t structure
 * \param[in]  cm             pointer to a cs_cell_mesh_t structure
 * \param[in]  pot            array of values of the potential (all the mesh)
 * \param[in, out] cb         auxiliary structure dedicated to diffusion
 * \param[in, out] vf_flux    array of values to set (size: n_vc)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_vbyf_flux(short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux);

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
cs_cdo_diffusion_wbs_face_flux(const cs_face_mesh_t      *fm,
                               const cs_real_t            pty_tens[3][3],
                               const double              *p_v,
                               const double               p_f,
                               const double               p_c,
                               cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal diffusive flux for a face assuming only the
 *          knowledge of the potential at faces and cell.
 *          CO+ST algorithm is used for reconstructing the normal flux from
 *          the degrees of freedom.
 *
 * \param[in]      f       face id in the cell mesh
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     array of values of the potential (all the mesh)
 * \param[in, out] cb      auxiliary structure dedicated to diffusion
 * \param[out]     flux    pointer to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_cost_flux(short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing a piecewise constant gradient from the degrees of
 *          freedom.
 *
 * \param[in]      f            face id in the cell mesh
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_tensor  property tensor times the face normal
 * \param[in]      pot_values   array of values of the potential (all the mesh)
 * \param[in, out] fluxes       values of the fluxes related to each vertex
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_p0_face_flux(const short int           f,
                                const cs_cell_mesh_t     *cm,
                                const cs_real_3_t        *diff_tensor,
                                const cs_real_t          *pot_values,
                                cs_real_t                *fluxes);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_DIFFUSION_H__ */
