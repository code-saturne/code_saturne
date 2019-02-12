#ifndef __CS_CDO_ADVECTION_H__
#define __CS_CDO_ADVECTION_H__

/*============================================================================
 * Build discrete convection operators for CDO schemes
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
#include "cs_property.h"

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
 * \brief   Define the local convection operator in CDO-Fb schemes
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes      array of advctive fluxes on primal faces
 * \param[in, out] adv         pointer to a cs_sdm_t structure to update
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_advection_t)(const cs_cell_mesh_t      *cm,
                       const cs_real_t            fluxes[],
                       cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_advection_bc_t)(const cs_equation_param_t   *eqp,
                          const cs_cell_mesh_t        *cm,
                          cs_cell_builder_t           *cb,
                          cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_advection_t)(const cs_equation_param_t   *eqp,
                       const cs_cell_mesh_t        *cm,
                       cs_real_t                    t_eval,
                       cs_face_mesh_t              *fm,
                       cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the advection operator in CDO
 *          vertex-based (or vertex+cell-based) schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] b       pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdovb_advection_bc_t)(const cs_cell_mesh_t       *cm,
                          const cs_equation_param_t  *eqp,
                          cs_real_t                   t_eval,
                          cs_face_mesh_t             *fm,
                          cs_cell_builder_t          *cb,
                          cs_cell_sys_t              *csys);

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the stabilization coefficient used in CIP scheme
 *
 * \param[in]  new_value     value of the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_set_cip_coef(double     new_value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the value of the stabilization coefficient used in CIP scheme
 *
 * \return   the value the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_advection_get_cip_coef(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the cellwise advection operator for CDO-Fb schemes
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval      time at which one evaluates the advection field
 * \param[in]      build_func  pointer to the function building the system
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_build(const cs_equation_param_t   *eqp,
                         const cs_cell_mesh_t        *cm,
                         cs_real_t                    t_eval,
                         cs_cdofb_advection_t        *build_func,
                         cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (without diffusion)
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc(const cs_equation_param_t   *eqp,
                       const cs_cell_mesh_t        *cm,
                       cs_cell_builder_t           *cb,
                       cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (with a diffusion term activated)
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc_wdi(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_cell_builder_t           *cb,
                           cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (without diffusion). Vector-valued case.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc_v(const cs_equation_param_t   *eqp,
                         const cs_cell_mesh_t        *cm,
                         cs_cell_builder_t           *cb,
                         cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (with a diffusion term activated). Vector-valued
 *          case.
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc_wdi_v(const cs_equation_param_t   *eqp,
                             const cs_cell_mesh_t        *cm,
                             cs_cell_builder_t           *cb,
                             cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the conservative formulation
 *         - upwind scheme
 *         - no diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] adv       pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwcsv(const cs_cell_mesh_t      *cm,
                           const cs_real_t            fluxes[],
                           cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the conservative formulation
 *         - upwind scheme
 *         - diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] adv       pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwcsv_di(const cs_cell_mesh_t      *cm,
                              const cs_real_t            fluxes[],
                              cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the non-conservative formulation
 *         - upwind scheme
 *         - no diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] cb        pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwnoc(const cs_cell_mesh_t      *cm,
                           const cs_real_t            fluxes[],
                           cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the non-conservative formulation
 *         - upwind scheme
 *         - diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] cb        pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwnoc_di(const cs_cell_mesh_t      *cm,
                              const cs_real_t            fluxes[],
                              cs_sdm_t                  *adv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when diffusion is activated and an upwind
 *          scheme and a conservative formulation is used
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv_di(const cs_equation_param_t   *eqp,
                              const cs_cell_mesh_t        *cm,
                              cs_real_t                    t_eval,
                              cs_face_mesh_t              *fm,
                              cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion and an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cencsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a mixed centered/upwind scheme with
 *          a conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_mcucsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when diffusion is activated and an upwind
 *          scheme and a conservative formulation is used
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc_di(const cs_equation_param_t   *eqp,
                              const cs_cell_mesh_t        *cm,
                              cs_real_t                    t_eval,
                              cs_face_mesh_t              *fm,
                              cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion when an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a non-conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cennoc(const cs_equation_param_t    *eqp,
                           const cs_cell_mesh_t         *cm,
                           cs_real_t                     t_eval,
                           cs_face_mesh_t               *fm,
                           cs_cell_builder_t            *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme when the advection field is cellwise
 *          constant
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_cw_cst(const cs_equation_param_t   *eqp,
                            const cs_cell_mesh_t        *cm,
                            cs_real_t                    t_eval,
                            cs_face_mesh_t              *fm,
                            cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb(const cs_equation_param_t   *eqp,
                     const cs_cell_mesh_t        *cm,
                     cs_real_t                    t_eval,
                     cs_face_mesh_t              *fm,
                     cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_bc(const cs_cell_mesh_t       *cm,
                       const cs_equation_param_t  *eqp,
                       cs_real_t                   t_eval,
                       cs_face_mesh_t             *fm,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator with CDO
 *          V+C schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_bc(const cs_cell_mesh_t        *cm,
                        const cs_equation_param_t   *eqp,
                        cs_real_t                    t_eval,
                        cs_face_mesh_t              *fm,
                        cs_cell_builder_t           *cb,
                        cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value of the upwinding coefficient in each cell knowing
 *          the related Peclet number
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in]      scheme    type of scheme used for the advection term
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_cell_upwind_coef(const cs_cdo_quantities_t    *cdoq,
                                  cs_param_advection_scheme_t   scheme,
                                  cs_real_t                     coefval[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_ADVECTION_H__ */
